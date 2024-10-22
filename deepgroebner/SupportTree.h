
#ifndef SUPPORT_TREE_H
#define SUPPORT_TREE_H

#include <vector>
#include <cstdint>

class SupportNode;

typedef std::vector<const int*> MoveList;

class SupportTree
{
public:
    SupportTree(int dimension);
    ~SupportTree();

    void add(const int32_t* move);
    void remove(const int32_t* move);
    void clear();

    const int32_t* reducable(
                    const int* move,
                    const int* skip_move = 0) const;
    const int32_t* reducable_negative(
                    const int* move,
                    const int* skip_move = 0) const;
    //void reducable( const Binomial& b,
    //                BinomialList& reducers) const;

    void print() const;

protected:
    const int32_t* reducable(
                    const int32_t* move,
                    const int32_t* skip_move,
                    const SupportNode* node) const;
    const int32_t* reducable_negative(
                    const int32_t* move,
                    const int32_t* skip_move,
                    const SupportNode* node) const;
    //void reducable( const Binomial& b,
    //                BinomialList& reducers,
    //                const SupportNode* node) const;

    inline bool reduces(const int32_t* move1, const int32_t* move2) const;
    inline bool reduces_negative(const int32_t* move1, const int32_t* move2) const;
		inline bool equals(const int32_t* move1, const int32_t* move2) const;

    void print(const SupportNode* node, std::vector<int> &labels) const;
		inline void print_move(const int32_t* move) const;

    SupportNode* root;
		int dim;
};

#endif
