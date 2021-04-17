import numpy as np
import tensorflow as tf

from deepgroebner.networks import ParallelEmbeddingLayer, ParallelDecidingLayer, Score, TransformerLayer, SelfAttentionLayer

class SelfAttentionLayer_Score_Q(tf.keras.layers.Layer):
    """A multi head self attention layer.

    Adapted from https://www.tensorflow.org/tutorials/text/transformer.

    Parameters
    ----------
    dim : int
        Positive integer dimension.
    n_heads : int, optional
        Positive integer number of heads (must divide `dim`).

    """

    def __init__(self, dim, softmax=True, n_heads=1):
        super(SelfAttentionLayer_Score_Q, self).__init__()
        assert dim % n_heads == 0, "number of heads must divide dimension"
        self.dim = dim
        self.n_heads = n_heads
        self.depth = dim // n_heads
        self.Wq = tf.keras.layers.Dense(dim)
        self.Wk = tf.keras.layers.Dense(dim)
        self.Wv = tf.keras.layers.Dense(dim)
        self.dense = tf.keras.layers.Dense(dim)
        self.attention_function=tf.nn.softmax
        if not softmax:
            self.score_attention = tf.keras.activations.sigmoid
        self.softmax = softmax
        self.supports_masking = True

        self.qval_learner = tf.keras.layers.Dense(dim)

        layers = [tf.keras.layers.Dense(128, activation = 'relu'),
                        tf.keras.layers.Dense(64, activation = 'relu'),
                        tf.keras.layers.Dense(32, activation = 'relu'),
                        tf.keras.layers.Dense(12, activation = 'relu'),
                        tf.keras.layers.Dense(1, activation = 'relu')]
        self.scorer = tf.keras.Sequential(layers=layers)

    def call(self, batch, mask=None):
        """Return the processed batch.

        Parameters
        ----------
        batch : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim, dim)
            Input batch with attached mask indicating valid rows.

        Returns
        -------
        output : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim, dim)
            Processed batch with mask passed through.

        """
        batch_size = tf.shape(batch)[0]
        Q = self.split_heads(self.Wq(batch), batch_size)
        K = self.split_heads(self.Wk(batch), batch_size)
        V = self.split_heads(self.Wv(batch), batch_size)
        Q_val = self.split_heads(self.get_qval(batch_size), batch_size)

        mask = mask[:, tf.newaxis, tf.newaxis, :]

        X = self.finish_attn(Q,K,V,batch_size,self.attention_function,mask=mask)
        if self.softmax:
            Y = self.finish_attn(Q_val,K,V,batch_size,self.attention_function,mask=mask)
        else:
            Y = self.finish_attn(Q_val,K,V,batch_size,self.score_attention, mask = mask)

        output = self.dense(X)
        score = self.scorer(Y)

        return output, score

    def get_qval(self, batch_size):
        return self.qval_learner(tf.ones([batch_size, 1, 1]))
    
    def finish_attn(self, Q, K, V, batch_size, attention_function, mask = None):
        X, attn_weights = self.scaled_dot_product_attention(Q, K, V, attention_function, mask=mask)
        X = tf.transpose(X, perm=[0, 2, 1, 3])
        X = tf.reshape(X, (batch_size, -1, self.dim))
        return X

    def split_heads(self, batch, batch_size):
        """Return batch reshaped for multihead attention."""
        X = tf.reshape(batch, (batch_size, -1, self.n_heads, self.depth))
        return tf.transpose(X, perm=[0, 2, 1, 3])

    def scaled_dot_product_attention(self, Q, K, V, attention_function, mask=None):
        """Return calculated vectors and attention weights.

        Parameters
        ----------
        Q : `Tensor` of type `tf.float32' and shape (..., dq, d1)
            Tensor of queries as rows.
        K : `Tensor` of type `tf.float32` and shape (..., dkv, d1)
            Tensor of keys as rows.
        V : `Tensor` of type `tf.float32` and shape (..., dkv, d2)
            Tensor of values as rows.
        mask : `Tensor of type `tf.bool' and shape (..., 1, dkv)
            The mask representing valid key/value rows.

        Returns
        -------
        output : `Tensor` of type `tf.float32` and shape (..., dq, d2)
            Processed batch of Q, K, V.
        attention_weights : `Tensor` of type `tf.float32` and shape (..., dq, dkv)
            Attention weights from intermediate step.

        """
        QK = tf.matmul(Q, K, transpose_b=True)
        d = tf.cast(tf.shape(K)[-1], tf.float32)
        attention_logits = QK / tf.math.sqrt(d)
        if mask is not None:
            attention_logits += tf.cast(~mask, tf.float32) * -1e9
        attention_weights = attention_function(attention_logits)
        output = tf.matmul(attention_weights, V)
        return output, attention_weights

class TransformerLayer_Score_Q(tf.keras.layers.Layer):
    """A transformer encoder layer.

    Parameters
    ----------
    dim : int
        Positive integer dimension of the attention layer and output.
    hidden_dim : int
        Positive integer dimension of the feed forward hidden layer.
    n_heads : int, optional
        Positive integer number of heads in attention layer (must divide `dim`).
    dropout : float, optional
        Dropout rate.

    """

    def __init__(self, dim, hidden_dim, softmax = True, n_heads=1, dropout=0.1):
        super(TransformerLayer_Score_Q, self).__init__()
        self.attention = SelfAttentionLayer_Score_Q(dim, softmax=softmax, n_heads=n_heads)
        self.dense1 = tf.keras.layers.Dense(hidden_dim, activation='relu')
        self.dense2 = tf.keras.layers.Dense(dim)
        self.dropout1 = tf.keras.layers.Dropout(dropout)
        self.dropout2 = tf.keras.layers.Dropout(dropout)
        self.layer_norm1 = tf.keras.layers.LayerNormalization(epsilon=1e-6)
        self.layer_norm2 = tf.keras.layers.LayerNormalization(epsilon=1e-6)
        self.supports_masking = True

    def call(self, batch, mask=None, training=False):
        """Return the processed batch.

        Parameters
        ----------
        batch : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim, dim)
            Input batch with attached mask indicating valid rows.

        Returns
        -------
        output : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim, dim)
            Processed batch with mask passed through.

        """
        X1, score= self.attention(batch, mask=mask)
        X1 = self.dropout1(X1, training=training)
        X1 = self.layer_norm1(batch + X1)
        X2 = self.dense2(self.dense1(X1))
        X2 = self.dropout2(X2, training=training)
        output = self.layer_norm2(X1 + X2)
        return output, score

class TransformerPMLP_Score_Q(tf.keras.Model):
    """A parallel multilayer perceptron network with a transformer layer.

    This model expects an input with shape (batch_dim, padded_dim, feature_dim), where
    entries are non-negative integers and padding is by -1. It returns a tensor
    of shape (batch_dim, padded_dim) where each batch is a softmaxed distribution over the rows
    with zero probability on any padded row.

    Parameters
    ----------
    dim : int
        Positive integer dimension of the transformer attention layer.
    hidden_dim : int
        Positive integer dimension of the transformer hidden feedforward layer.
    activation : {'relu', 'selu', 'elu', 'tanh', 'sigmoid'}, optional
        Activation for the embedding.
    final_activation : {'log_softmax', 'softmax'}, optional
        Activation for the final output layer.

    """

    def __init__(self, dim, hidden_dim, softmax = True, activation='relu', final_activation='log_softmax'):
        super(TransformerPMLP_Score_Q, self).__init__()
        self.embedding = ParallelEmbeddingLayer(dim, [], final_activation=activation)
        self.attn = TransformerLayer_Score_Q(dim, hidden_dim, softmax, n_heads=4)
        self.deciding = ParallelDecidingLayer([], final_activation=final_activation)

    def call(self, batch):
        X = self.embedding(batch)
        X, score = self.attn(X)
        X = self.deciding(X)
        return X, tf.squeeze(score, axis = -1)

class SelfAttentionLayer_DVal(SelfAttentionLayer):
    """A multi head self attention layer.

    Adapted from https://www.tensorflow.org/tutorials/text/transformer.

    Parameters
    ----------
    dim : int
        Positive integer dimension.
    n_heads : int, optional
        Positive integer number of heads (must divide `dim`).

    """

    def __init__(self, dim, softmax=True, n_heads=1):
        super().__init__(dim, n_heads=n_heads)
        if softmax:
            self.attention_function= tf.nn.softmax
        else:
            self.attention_function = tf.keras.activations.sigmoid
        self.qval_learner = tf.keras.layers.Dense(dim)
        layers = [tf.keras.layers.Dense(128, activation = 'relu'),
                        tf.keras.layers.Dense(64, activation = 'relu'),
                        tf.keras.layers.Dense(32, activation = 'relu'),
                        tf.keras.layers.Dense(12, activation = 'relu'),
                        tf.keras.layers.Dense(1, activation = 'relu')]
        self.scorer = tf.keras.Sequential(layers=layers)

    def call(self, batch, mask=None):
        """Return the processed batch.

        Parameters
        ----------
        batch : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim, dim)
            Input batch with attached mask indicating valid rows.

        Returns
        -------
        output : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim, dim)
            Processed batch with mask passed through.

        """
        batch_size = tf.shape(batch)[0]
        K = super().split_heads(self.Wk(batch), batch_size)
        V = super().split_heads(self.Wv(batch), batch_size)
        Q_val = super().split_heads(self.get_qval(batch_size), batch_size)

        mask = mask[:, tf.newaxis, tf.newaxis, :]

        Y = self.finish_attn(Q_val,K,V,batch_size,mask=mask)
        score = self.scorer(Y)
        return tf.math.add(score, tf.constant(1, dtype = tf.float32))

    def get_qval(self, batch_size):
        return self.qval_learner(tf.ones([batch_size, 1, 1]))
    
    def finish_attn(self, Q, K, V, batch_size, mask = None):
        X, attn_weights = super().scaled_dot_product_attention(Q, K, V, mask=mask)
        X = tf.transpose(X, perm=[0, 2, 1, 3])
        X = tf.reshape(X, (batch_size, -1, self.dim))
        return X

class TransformerLayer_DVal(TransformerLayer):
    """A transformer encoder layer.

    Parameters
    ----------
    dim : int
        Positive integer dimension of the attention layer and output.
    hidden_dim : int
        Positive integer dimension of the feed forward hidden layer.
    n_heads : int, optional
        Positive integer number of heads in attention layer (must divide `dim`).
    dropout : float, optional
        Dropout rate.

    """

    def __init__(self, dim, hidden_dim, n_heads=1, dropout=0.1):
        super().__init__(dim, hidden_dim, n_heads, dropout)
        self.attention = SelfAttentionLayer_Score_Q(dim, softmax=True, n_heads=n_heads)
        self.attention_v2 = SelfAttentionLayer_DVal(dim, softmax = False, n_heads=n_heads)


    def call(self, batch, mask=None, training=False):
        """Return the processed batch.

        Parameters
        ----------
        batch : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim, dim)
            Input batch with attached mask indicating valid rows.

        Returns
        -------
        output : `Tensor` of type `tf.float32` and shape (batch_dim, padded_dim, dim)
            Processed batch with mask passed through.

        """
        X1, score_softmax= self.attention(batch, mask=mask)
        score_sigmoid= self.attention_v2(batch, mask = mask)
        X1 = self.dropout1(X1, training=training)
        X1 = self.layer_norm1(batch + X1)
        X2 = self.dense2(self.dense1(X1))
        X2 = self.dropout2(X2, training=training)
        output = self.layer_norm2(X1 + X2)
        return output, score_softmax, score_sigmoid

class Score_Decider(tf.keras.layers.Layer):

    def __init__(self, hidden_layers = [], last_activation = 'relu'):
        super(Score_Decider, self).__init__()
        layers = []
        for dim in hidden_layers:
            layers.append(tf.keras.layers.Dense(dim, activation = 'elu'))
        layers.append(tf.keras.layers.Dense(1, activation=last_activation))
        self.scorer = tf.keras.Sequential(layers=layers)
    
    def call(self, batch):
        X = self.scorer(batch)
        return X

class TransformerPMLP_DVal(tf.keras.Model):
    def __init__(self, dim, hidden_dim, hidden_layers = [], activation='relu', score_activation = 'relu', final_activation='log_softmax'):
        super(TransformerPMLP_DVal, self).__init__()
        self.embedding = ParallelEmbeddingLayer(dim, [], final_activation=activation)
        self.attn = TransformerLayer_DVal(dim, hidden_dim, n_heads=4)
        self.scorer = Score_Decider(hidden_layers, score_activation)
        self.deciding = ParallelDecidingLayer([], final_activation=final_activation)
        
    def call(self, batch):
        X = self.embedding(batch)
        X, score_softmax, score_sigmoid = self.attn(X)
        score = self.scorer(tf.concat([score_softmax, score_sigmoid], axis = len(score_softmax.shape)-1)) # think how to do this better
        X = self.deciding(X)
        return X, score

class TransformerPMLP_MHA_Q_Scorer(tf.keras.Model):
    def __init__(self, dim, hidden_dim, hidden_layers = [], activation='relu', score_activation = 'relu', final_activation='log_softmax'):
        super(TransformerPMLP_MHA_Q_Scorer, self).__init__()
        self.embedding = ParallelEmbeddingLayer(dim, [], final_activation=activation)
        self.attn = TransformerLayer_Score_Q(dim, hidden_dim, True, n_heads=4)
        self.scorer = Score(hidden_layers)
        self.score_blender = Score_Decider(hidden_layers)
        self.deciding = ParallelDecidingLayer([], final_activation=final_activation)
    
    def call(self, batch):
        X = self.embedding(batch)
        X, score_q = self.attn(X)
        score_mha = self.scorer(X)
        score = self.score_blender(tf.concat([score_q, score_mha], axis = len(score_mha.shape)-1))
        X = self.deciding(X)
        return X, score
