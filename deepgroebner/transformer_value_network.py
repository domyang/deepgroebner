"""
Transformer that is just a value network
"""

from deepgroebner.networks import ParallelEmbeddingLayer, ParallelDecidingLayer

import numpy as np
import tensorflow as tf

class SelfAttentionLayer(tf.keras.layers.Layer):

    """A multi head self attention layer.

    Adapted from https://www.tensorflow.org/tutorials/text/transformer.

    Parameters
    ----------
    dim : int
        Positive integer dimension.
    n_heads : int, optional
        Positive integer number of heads (must divide `dim`).

    """

    def __init__(self, dim, softmax = True, n_heads=1):
        super(SelfAttentionLayer, self).__init__()
        assert dim % n_heads == 0, "number of heads must divide dimension"
        self.dim = dim
        self.n_heads = n_heads
        self.depth = dim // n_heads
        #self.Wq = tf.keras.layers.Dense(dim)
        self.Wk = tf.keras.layers.Dense(dim)
        self.Wv = tf.keras.layers.Dense(dim)
        self.dense = tf.keras.layers.Dense(dim)
        self.Value_Q = self.set_q(dim)
        if softmax:
            self.attention_function = tf.nn.softmax
        else:
            self.attention_function = tf.nn.sigmoid
        self.supports_masking = True

    def set_q(self, dim):
        temp = tf.random.uniform(shape=(1, 1, dim), dtype = tf.float32, seed = 0)
        return tf.Variable(temp, trainable = True)
   
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
        Q = self.split_heads(tf.tile(self.Value_Q, [batch_size, 1, 1]), batch_size)
        K = self.split_heads(self.Wk(batch), batch_size)
        V = self.split_heads(self.Wv(batch), batch_size)
        mask = mask[:, tf.newaxis, tf.newaxis, :]
        X, attn_weights = self.scaled_dot_product_attention(Q, K, V, mask=mask) 
        X = tf.transpose(X, perm=[0, 2, 1, 3])
        output = tf.reshape(X, (batch_size, -1, self.dim)) # (batch_dim, 1, dim)
        return output

    def split_heads(self, batch, batch_size):
        """Return batch reshaped for multihead attention."""
        X = tf.reshape(batch, (batch_size, -1, self.n_heads, self.depth))
        return tf.transpose(X, perm=[0, 2, 1, 3])

    def scaled_dot_product_attention(self, Q, K, V, mask=None):
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
        attention_weights = tf.math.add(self.attention_function(attention_logits), tf.constant(1, dtype=tf.float32))
        output = tf.matmul(attention_weights, V)
        return output, attention_weights

class Value_Function(tf.keras.layers.Layer):
    """
    Value Function to predict expected discounted rewards based on output of multi-headed-attention layer
    """

    def __init__(self, hidden_layers:list):
        super(Value_Function, self).__init__()
        self.value_function = tf.keras.Sequential()
        for layer in hidden_layers:
            self.value_function.add(tf.keras.layers.Dense(layer, activation='relu'))
        self.value_function.add(tf.keras.layers.Dense(1, activation='relu'))

    def call(self, batch):
        return tf.math.negative(self.value_function(batch))

class TransformerValueModel(tf.keras.models.Model):
    """
    Combination of the self attention layer of a transformer and value function
    """
    
    def __init__(self, hidden_layers, dim, softmax, n_heads = 4):
        super(TransformerValueModel, self).__init__()
        self.embedding = ParallelEmbeddingLayer(dim, [])
        self.value_model = Value_Function(hidden_layers)
        self.mha = SelfAttentionLayer(dim, softmax, n_heads)

    def call(self, batch):
        X = self.embedding(batch)
        X = self.mha(X)
        value = self.value_model(X)
        return tf.squeeze(value, axis = -1)
