"""Neural networks for agents."""

import numpy as np
import tensorflow as tf


class MultilayerPerceptron(tf.keras.Model):
    """A basic multilayer perceptron network.

    This network is used for the policy and value models when training on CartPole-v0,
    CartPole-v1, and LunarLander-v2.

    Parameters
    ----------
    output_dim : int
        Positive integer output dimension of the network.
    hidden_layers : list
        List of positive integer hidden layer dimensions.
    activation : {'relu', 'selu', 'elu', 'tanh', 'sigmoid'}, optional
        Activation used for the hidden layers.
    final_activation : {'log_softmax', 'softmax', 'linear', 'exponential'}, optional
        Activation used for the final output layer.

    Examples
    --------
    >>> mlp = MultilayerPerceptron(2, [128])
    >>> states = tf.random.uniform((64, 4))
    >>> logprobs = mlp(states)
    >>> logprobs.shape
    TensorShape([64, 2])
    >>> actions = tf.random.categorical(logprobs, 1)
    >>> actions.shape
    TensorShape([64, 1])

    """

    def __init__(self, output_dim, hidden_layers, activation='relu', final_activation='log_softmax'):
        super(MultilayerPerceptron, self).__init__()
        final_activation = tf.nn.log_softmax if final_activation == 'log_softmax' else final_activation
        self.hidden_layers = [tf.keras.layers.Dense(u, activation=activation) for u in hidden_layers]
        self.final_layer = tf.keras.layers.Dense(output_dim, activation=final_activation)

    def call(self, X):
        for layer in self.hidden_layers:
            X = layer(X)
        return self.final_layer(X)


class ParallelEmbeddingLayer(tf.keras.layers.Layer):
    """A layer for computing a nonlinear embedding of non-negative integer feature vectors.

    This layer is used with the LeadMonomialsWrapper to embed the exponent vectors of pairs
    into feature vectors. Each vector is embedded independently by a single learned multilayer
    perceptron. A mask is generated and attached to the output based on padding by -1.

    Parameters
    ----------
    embed_dim : int
        Positive integer output dimension of the embedding.
    hidden_layers : list
        List of positive integer hidden layer dimensions.
    activation : {'relu', 'selu', 'elu', 'tanh', 'sigmoid'}, optional
        Activation used for the hidden layers.
    final_activation : {'relu', 'linear', 'exponential'}, optional
        Activation used for the final output embedding layer.

    """

    def __init__(self, embed_dim, hidden_layers, activation='relu', final_activation='relu'):
        super(ParallelEmbeddingLayer, self).__init__()
        self.hidden_layers = [tf.keras.layers.Dense(u, activation=activation) for u in hidden_layers]
        self.final_layer = tf.keras.layers.Dense(embed_dim, activation=final_activation)

    def call(self, batch):
        """Return the embedding for this batch.

        Parameters
        ----------
        batch : `Tensor` of type `tf.int32` and shape (batch_dim, padded_dim, input_dim)
            Input batch, with padded rows indicated by -1 and all other values non-negative.

        Returns
        -------
        output : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim, embed_dim)
            Embedding of the input batch with attached mask indicating valid rows.

        """
        X = tf.cast(batch, tf.float32)
        for layer in self.hidden_layers:
            X = layer(X)
        output = self.final_layer(X)
        return output

    def compute_mask(self, batch, mask=None):
        return tf.math.not_equal(batch[:, :, -1], -1)


class RecurrentEmbeddingLayer(tf.keras.layers.Layer):
    """A layer for computing a nonlinear embedding of non-negative integer feature vectors.

    This layer is used with the LeadMonomialsWrapper to embed the exponent vectors of pairs
    into feature vectors. An RNN is used so vector embeddings can depend on other vectors.
    A mask is generated and attached to the output based on padding by -1.

    Parameters
    ----------
    embed_dim : int
        Positive integer output dimension of the embedding.
    hidden_layers : list
        List of positive integer hidden layer dimensions.
    cell : {'gru', 'lstm'}
        Type of recurrent cell.

    """

    def __init__(self, embed_dim, hidden_layers, cell='gru'):
        super(RecurrentEmbeddingLayer, self).__init__()
        self.cell = cell
        cell_fn = tf.keras.layers.GRU if cell == 'gru' else tf.keras.layers.LSTM
        self.hidden_layers = [cell_fn(u, return_sequences=True) for u in hidden_layers]
        self.final_layer = cell_fn(embed_dim, return_sequences=True, return_state=True)

    def call(self, batch):
        """Return embedding and final hidden states for this batch.

        Parameters
        ----------
        batch : `Tensor` of type tf.int32` and shape (batch_dim, padded_dim, feature_dim)
            Input batch, with padded rows indicated by -1 and all other values non-negative.

        Returns
        -------
        output : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim, embed_dim)
            Embedding of the input batch with mask indicating valid rows.
        state : `Tensor` or pair of `Tensor`s of type `tf.float32` and shape (batch_dim, embed_dim)

        """
        mask = self.compute_mask(batch)
        X = tf.cast(batch, tf.float32)
        for layer in self.hidden_layers:
            X, *state = layer(X, mask=mask)
        output, *state = self.final_layer(X, mask=mask)
        return (output, *state)

    def compute_mask(self, batch, mask=None):
        return tf.math.not_equal(batch[:, :, -1], -1)


class DenseProcessingLayer(tf.keras.layers.Layer):
    """A simple processing stack that applies dense layers.

    Parameters
    ----------
    output_dim : int
       Positive integer output dimension of this layer.
    hidden_layers : list
       List of positive integer hidden layer dimensions
    activation : {'relu', 'selu', 'elu', 'tanh', 'sigmoid'}, optional
        Activation used for the hidden layers and output.    

    """

    def __init__(self, output_dim, hidden_layers, activation='relu'):
        super(DenseProcessingLayer, self).__init__()
        self.hidden_layers = [tf.keras.layers.Dense(u, activation=activation) for u in hidden_layers]
        self.final_layer = tf.keras.layers.Dense(output_dim, activation=activation)
        self.supports_masking = True

    def call(self, batch):
        """Return the processed output for this batch.

        Parameters
        ----------
        batch : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim, input_dim)
            Input batch with attached mask indicating valid rows.

        Returns
        -------
        output : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim, output_dim)
            Processed batch with attached mask passed through.

        """
        X = batch
        for layer in self.hidden_layers:
            X = layer(X)
        output = self.final_layer(X)
        return output


class SelfAttentionLayer(tf.keras.layers.Layer):
    """A multi head self attention layer.

    Adapted from https://www.tensorflow.org/tutorials/text/transformer.

    """

    def __init__(self, dim, n_heads):
        super(MultiHeadSelfAttentionLayer, self).__init__()
        assert dim % n_heads == 0
        self.dim = dim
        self.n_heads = n_heads
        self.depth = dim // n_heads
        self.Wq = tf.keras.layers.Dense(dim)
        self.Wk = tf.keras.layers.Dense(dim)
        self.Wv = tf.keras.layers.Dense(dim)
        self.dense = tf.keras.layers.Dense(dim)

    def call(self, X, mask=None):
        batch_size = tf.shape(X)[0]
        Q = self.split_heads(self.Wq(X), batch_size)
        K = self.split_heads(self.Wk(X), batch_size)
        V = self.split_heads(self.Wv(X), batch_size)
        mask = mask[:, tf.newaxis, tf.newaxis, :]
        scaled_attention, attention_weights = scaled_dot_product_attention(Q, K, V, mask=mask)
        scaled_attention = tf.transpose(scaled_attention, perm=[0, 2, 1, 3])
        concat_attention = tf.reshape(scaled_attention, (batch_size, -1, self.dim))
        X = self.dense(concat_attention)
        return X, attention_weights

    def split_heads(self, X, batch_size):
        X = tf.reshape(X, (batch_size, -1, self.n_heads, self.depth))
        return tf.transpose(X, perm=[0, 2, 1, 3])

    def scaled_dot_product_attention(Q, K, V, mask=None):
        """Return calculated vector and attention weights.

        Parameters
        ----------
        Q : tensor of shape (..., dq, d)
            Tensor of queries as rows.
        K : tensor of shape (..., dk, d)
            Tensor of keys as rows.
        V : tensor of shape (..., dk, dv)
            Tensor of values as rows.
        mask : boolean tensor of shape broadcastable to (..., dq, dk)
            The mask representing valid rows.

        Returns
        -------
        X : tensor
        attention_weights : tensor

        """
        QK = tf.matmul(Q, K, transpose_b=True)
        d = tf.cast(tf.shape(K)[-1], tf.float32)
        attention_logits = QK / tf.math.sqrt(d)
        if mask is not None:
            attention_logits += (mask * -1e9)
        attention_weights = tf.nn.softmax(attention_logits)
        X = tf.matmul(attention_weights, V)
        return X, attention_weights


class TransformerLayer(tf.keras.layers.Layer):

    def __init__(self, num_heads, input_dim, feed_forward_hidden_size, training: bool, rate = .1):
        '''
        Constructor.

        @Params:
            num_heads - Number of self attention layers NOTE: num_heads must divide input_dim
            input_dim - Dimension of input 
            feed_forward_hidden_size - size of feed forward network at the end of the encoder
            training - encorporate drop out 
            rate - rate of dropout
        '''
        super(EncoderLayer, self).__init__()
        self.attention = MultiHeadSelfAttentionLayer(input_dim, num_heads)
        self.first_lt = tf.keras.layers.Dense(feed_forward_hidden_size, activation = 'relu')
        self.second_lt = tf.keras.layers.Dense(input_dim)

        self.layer_norm_mha = tf.keras.layers.LayerNormalization(epsilon=1e-6)
        self.layer_norm_ff = tf.keras.layers.LayerNormalization(epsilon=1e-6)

        self.dropout1 = tf.keras.layers.Dropout(rate)
        self.dropout2 = tf.keras.layers.Dropout(rate)

        self.training = training

    def call(self, input_set):

        input_set = tf.cast(input_set, tf.float32)

        # Multiheaded attention
        attention_output = self.attention(input_set, input_set, input_set)
        attention_output = self.dropout1(attention_output, training= self.training) # Only use drop out when training
        att_norm = self.layer_norm_mha(input_set + attention_output) # Residual connection and norm

        # Feed forward network stage
        ff = self.second_lt(self.first_lt(att_norm))
        ff1 = self.dropout1(ff, training = self.training)
        encoder_output = self.layer_norm_ff(input_set + ff1)
        return encoder_output


class ParallelDecidingLayer(tf.keras.layers.Layer):
    """A layer for computing probability distributions over arbitrary numbers of rows.

    This layer is used following an embedding and processing of the state of a LeadMonomialsWrapper
    to produce the policy probabilities for each available action. The layer learns a single function
    implemented with a multilayer perceptron that computes scores independently for each row. These
    scores are softmaxed to produce the probabilites.

    Parameters
    ----------
    hidden_layers : list
        List of positive integer hidden layer dimensions.
    activation : {'relu', 'selu', 'elu', 'tanh', 'sigmoid'}, optional
        Activation for the hidden layers.
    final_activation : {'log_softmax', 'softmax'}, optional
        Activation for the final output layer.

    """

    def __init__(self, hidden_layers, activation='relu', final_activation='log_softmax'):
        super(ParallelDecidingLayer, self).__init__()
        self.hidden_layers = [tf.keras.layers.Dense(u, activation=activation) for u in hidden_layers]
        self.final_layer = tf.keras.layers.Dense(1, activation='linear')
        self.final_activation = tf.nn.log_softmax if final_activation == 'log_softmax' else tf.nn.softmax

    def call(self, batch, mask=None):
        """Return probability distributions over rows of batch.

        Parameters
        ----------
        batch : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim, feature_dim)
            Batch of feature vectors with attached mask indicating valid rows.

        Returns
        -------
        output : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim)
            Softmaxed probability distributions over valid rows.

        """
        X = batch
        for layer in self.hidden_layers:
            X = layer(X)
        if mask is not None:
            X = tf.squeeze(self.final_layer(X), axis=-1) - 1e9 * tf.cast(~mask, tf.float32)
        else:
            X = tf.squeeze(self.final_layer(X), axis=-1)
        output = self.final_activation(X)
        return output


class PointerDecidingLayer(tf.keras.layers.Layer):
    """Pointer network attention mechanism. This will also handle the one decode step
    
    input shape is (batch_dim, padded_dim, feature_dim)
    output shape is (batch_dim, padded_dim)

        Attribute -
            decoder_lstm: Only used for starting token which is a fixed random vector
            encoder_weight: Matrix to project encoder output during attention
            decode_weight: Matrix to project decode output (the start token)
            v: Used to project encode and decode into a 1d number
            tanh: activation function for attention
            input_size: size of input
            hidden_size: number of hidden nodes
    """

    def __init__(self, input_dim, hidden_layer, layer_type = 'lstm', dot_product_attention = False, prob = 'norm'):
        super(pointer, self).__init__()
        if layer_type == 'lstm':
            self.decoder_layer = tf.keras.layers.LSTM(hidden_layer, return_sequences=True)
        elif layer_type == 'gru':
            self.decoder_layer = tf.keras.layers.GRU(hidden_layer, return_sequences=True) 
        if not dot_product_attention:
            self.encoder_weight = tf.keras.layers.Dense(hidden_layer)
            self.decode_weight = tf.keras.layers.Dense(hidden_layer)
            self.v = tf.keras.layers.Dense(1)
            self.tanh = tf.keras.activations.tanh
        if prob == 'log':
            self.softmax = tf.nn.log_softmax
        else:
            self.softmax = tf.nn.softmax
        self.dot_prod_attention = dot_product_attention
        self.input_size = input_dim
        self.hidden_size = hidden_layer

    def call(self, encoder_output, initial_states=None):
        '''
        Calculate attention.
        
        Params:
            encoder_output: output of the encoder block
                size: (batch, seq_len, input_dim)
            initial_states: hidden and cell states from the encoder block
       '''
        batch_size = encoder_output.shape[0]
        start_token = self.initialize_start_token(batch_size)
        lstm_decoder_output = self.decoder_layer(start_token, initial_state=initial_states) #(batch, 1, input_dim)
        if self.dot_prod_attention:
            attention_scores = tf.squeeze(tf.linalg.matmul(lstm_decoder_output, encoder_output, transpose_b = True), axis = 1)
            return self.softmax(attention_scores)
        else:
            similarity_score = np.zeros([batch_size, encoder_output.shape[1]])
            lstm_decoder_projection = self.decode_weight(lstm_decoder_output)
            encoder_project = self.encoder_weight(encoder_output)
            for batch in range(batch_size):
                for index in range(similarity_score.shape[1]):
                    #e_i = tf.expand_dims(encoder_output[batch][index], axis = 0)
                    similarity_score[batch][index] = self.v(self.tanh(encoder_project[batch][index] + lstm_decoder_projection[batch]))[0][0] # change this
            return self.softmax(tf.convert_to_tensor(similarity_score))

    def initialize_start_token(self, batch_size):
        '''
        Initialize start token

        Params:
            batch_size: size of batch
        '''
        np.random.seed(42)
        start_token = tf.convert_to_tensor(np.random.random([batch_size,1,self.input_size]).astype(np.float32)) # Change this
        np.random.seed()
        return start_token


class ParallelMultilayerPerceptron(tf.keras.Model):
    """A parallel multilayer perceptron network.

    This model expects an input with shape (batch_dim, padded_dim, feature_dim), where
    entries are non-negative integers and padding is by -1. It returns a tensor
    of shape (batch_dim, padded_dim) where each batch is a softmaxed distribution over the rows
    with zero probability on any padded row.

    Parameters
    ----------
    hidden_layers : list
        List of positive integer hidden layer dimensions.
    activation : {'relu', 'selu', 'elu', 'tanh', 'sigmoid'}, optional
        Activation for the hidden layers.
    final_activation : {'log_softmax', 'softmax'}, optional
        Activation for the final output layer.

    Examples
    --------
    >>> pmlp = ParallelMultilayerPerceptron([128])
    >>> states = tf.constant([
            [[ 0,  1],
             [ 3,  0],
             [-1, -1]],
            [[ 8,  5],
             [ 3,  3],
             [ 3,  5]],
            [[ 6,  7],
             [ 6,  8],
             [-1, -1]],
        ])
    >>> logprobs = pmlp(states)
    >>> logprobs.shape
    TensorShape([3, 3])
    >>> actions = tf.random.categorical(logprobs, 1)
    >>> actions.shape
    TensorShape([3, 1])

    """

    def __init__(self, hidden_layers, activation='relu', final_activation='log_softmax'):
        super(ParallelMultilayerPerceptron, self).__init__()
        self.embedding = ParallelEmbeddingLayer(hidden_layers[-1], hidden_layers[:-1],
                                                activation=activation, final_activation=activation)
        self.deciding = ParallelDecidingLayer([], final_activation=final_activation)

    def call(self, batch):
        X = self.embedding(batch)
        X = self.deciding(X)
        return X


class PointerNetwork(tf.keras.Model):
    """Recurrent embedding followed by pointer."""

    def __init__(self, input_dim, hidden_layer, input_layer = 'lstm', dot_prod_attention=False, prob = 'norm'):
        '''
        Params:
            input_dim: dimension of input
            hidden_layer: dimension of output layer
            input_layer: lstm or gru
            dot_prod_attention: dot product attention or traditional pointer network attention
        '''
        super(PointerNetwork, self).__init__()
        self.encoder = pnetEncoder(hidden_layer, layer_type = input_layer)
        self.point = pointer(input_dim, hidden_layer, layer_type=input_layer, dot_product_attention=dot_prod_attention, prob = prob)
        self.hidden_size = hidden_layer
        self.layer = input_layer
        self.attention_type = dot_prod_attention

    def __call__(self, input):
        '''
        '''
        if self.layer == 'lstm':
            seq_output, mem_state, carry_state = self.encoder(input)
            initial_states = [mem_state, carry_state]
        else:
            seq_output, mem_state = self.encoder(input)
            initial_states = mem_state
        prob_dist = self.point(seq_output, initial_states)

        return prob_dist


class CustomLSTM(tf.keras.layers.Layer):
    def __init__(self, hl_out_size):
        super(CustomLSTM, self).__init__()
        self.Whi = tf.keras.layers.Dense(hl_out_size, activation='sigmoid')
        self.Whf = tf.keras.layers.Dense(hl_out_size, activation='sigmoid')
        self.Whg = tf.keras.layers.Dense(hl_out_size, activation='tanh')
        self.Who = tf.keras.layers.Dense(hl_out_size, activation='sigmoid')
    def __call__(self, input, cell_state):
        i_t = self.Whi(input)
        f_t = self.Whf(input)
        g_t = self.Whg(input)
        o_t = self.Who(input)
        c_t = tf.math.multiply(f_t,cell_state) + tf.math.multiply(i_t, g_t)
        h_t = tf.math.multiply(o_t, tf.tanh(c_t))
        return h_t, c_t


class ProcessBlock(tf.keras.layers.Layer):
    def __init__(self, hidden_layer, num_step):
        '''
        Constructor for ProcessBlock.

        @Params:
            hidden_layer: Size of embeddings and LSTM
            num_step: Number of times I process the input
        '''
        super(ProcessBlock, self).__init__()
        self.embed = tf.keras.layers.Dense(hidden_layer)
        #self.convul = tf.keras.layers.Dense(hidden_layer)
        #self.process_block = tf.keras.layers.GRU(hidden_layer, return_state=True)
        self.process_block = CustomLSTM(hidden_layer)
        self.hidden_size = hidden_layer
        self.num_step = num_step
    def read_out(self, M, q, c,batch_size):
        '''
        Perform the attention with the memory vectors and pass by the LSTM

        @Params:
            M - Polynomial Embeddings
            q - Query vector (last hidden state)
            c - Last cell state
            batch_size - size of batch
        '''
        attention = tf.nn.softmax(tf.linalg.matmul(M, q, transpose_b = True))
        r_t = tf.linalg.matmul(attention, M, transpose_a = True) # sum of weight memory vector (by attention)
        q_star_t = tf.concat([q, r_t], axis = 2)
        mem_state, cell_state = self.process_block(q_star_t, c)
        return mem_state, cell_state
    def initHiddenState(self, batch_size):
        '''
        Initialize initial hidden states and cell state. 

        @Params:
            batch_size - size of batch
        '''
        np.random.seed(42)
        hidden_state = tf.convert_to_tensor(np.random.random([batch_size,1,self.hidden_size]).astype(np.float32)) # Change this
        cell_state = tf.convert_to_tensor(np.random.random([batch_size,1,self.hidden_size]).astype(np.float32)) # Change this
        np.random.seed()
        return hidden_state, cell_state
    def __call__(self, input_seq):
        '''
        Calculate embedding.

        @Params:
            input_seq - input sequence
        '''
        embeddings = self.embed(input_seq)
        initial_state, cell_state = self.initHiddenState(input_seq.shape[0])

        # Start processing
        for _ in range(self.num_step):
            initial_state, cell_state = self.read_out(embeddings, initial_state, cell_state, input_seq.shape[0])
        return initial_state, embeddings


class PBPointerNet(tf.keras.Model):
    def __init__(self, input_dim, hidden_layer, layer_type = 'lstm', dot_product_attention = False, prob = 'norm'):
        super(PBPointerNet, self).__init__()
        self.pointer = pointer(input_dim, hidden_layer, layer_type, dot_product_attention, prob)
        self.processBlock = ProcessBlock(hidden_layer, 6)

    def __call__(self, input):
        output, embeddings = self.processBlock(input)
        output = tf.squeeze(output, axis = 1)
        prob = self.pointer(embeddings, output)
        return prob


class PairsLeftBaseline:
    """A Buchberger value network that returns discounted pairs left."""

    def __init__(self, gam=0.99):
        self.gam = gam
        self.trainable_variables = []

    def predict(self, X, **kwargs):
        states, pairs, *_ = X.shape
        if self.gam == 1:
            fill_value = - pairs
        else:
            fill_value = - (1 - self.gam ** pairs) / (1 - self.gam)
        return np.full((states, 1), fill_value)

    def __call__(self, inputs):
        return self.predict(inputs)

    def save_weights(self, filename):
        pass

    def load_weights(self, filename):
        pass


class AgentBaseline:
    """A Buchberger value network that returns an agent's performance."""

    def __init__(self, agent, gam=0.99):
        self.agent = agent
        self.gam = gam
        self.trainable_variables = []

    def predict(self, env):
        env = env.copy()
        R = 0.0
        discount = 1.0
        state = (env.G, env.P) if hasattr(env, 'G') else env._matrix()
        done = False
        while not done:
            action = self.agent.act(state)
            state, reward, done, _ = env.step(action)
            R += reward * discount
            discount *= self.gam
        return R

    def __call__(self, inputs):
        return self.predict(inputs)

    def save_weights(self, filename):
        pass

    def load_weights(self, filename):
        pass
