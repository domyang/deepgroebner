
#include "SupportTree.h"

#include <vector>
#include <iostream>
#include <cassert>
#include <cstdint>


class SupportNode {
public:
    SupportNode() { moves = 0;}
    virtual ~SupportNode()
    {
        delete moves;
        for (int i = 0; i < (int) nodes.size(); ++i) { delete nodes[i].second; }
    }
    std::vector<std::pair<int,SupportNode*> > nodes;
    MoveList* moves;
    // Filter* filter;
};

SupportTree::SupportTree(int dimension)
{
    root = new SupportNode();
		dim = dimension;
}

SupportTree::~SupportTree()
{
    delete root;
}

void
SupportTree::add(const int32_t* move)
{
    SupportNode* current = root;
    for (int i = 0; i < dim; ++i) {

				// std::cout << "index " << i << std::endl;

				// If move has positive support at i
        if (move[i] > 0) {
            int j = 0;
						// If there is a node with label i in node list, go to this node
            while (j < (int) current->nodes.size() && current->nodes[j].first != i){
								++j;
						}

            if (j < (int) current->nodes.size()) {
								// std::cout << "Node " << j << std::endl;
                current = current->nodes[j].second;
            }
            else { // Otherwise, create a node with this label and go to this node
								// std::cout << "Creating new node" << std::endl;
                SupportNode* next = new SupportNode;
                current->nodes.push_back(std::pair<int,SupportNode*>(i,next));
                current = next;
            }
        }
    }
		// If this node has no moves in it, we create a MoveList and add this move to it
    if (!current->moves) {
				// std::cout << "Creating MoveList" << std::endl;
        current->moves = new MoveList;
        // current->filter = new Filter;
        // b.get_filter(*current->filter);
    }

		int32_t* move_copy = new int32_t[dim];
		for (int i = 0; i < dim; ++i)
				move_copy[i] = move[i];
    current->moves->push_back(move_copy);
		// std::cout << "Added!" << std::endl;
}

// Assumes point exists.
void SupportTree::remove(const int32_t* move)
{
    SupportNode* current = root;
    for (int i = 0; i < dim; ++i) {
				// If move supported at i
        if (move[i] > 0) {
            int j = 0;

						// Search 
            while (j < (int) current->nodes.size() && current->nodes[j].first != i)
								++j; 
						// We reached the end of the list and did not find
            if (j == (int) current->nodes.size()) {
                assert(false); // If we got to here, then we did not find the point.
            }
            else {
                current = current->nodes[j].second;
            }
        }
    }

		// Assert that this node has moves
    assert(current->moves); // If this assert is false, the point does not exist.
    for (MoveList::iterator iter = current->moves->begin();
                    iter != current->moves->end(); ++iter) {
				bool equal = true;
				for (int j = 0; j < dim; j++) {
						if ((*iter)[j] != move[j]) {
								equal = false;
								break;
						}
				}

        if (equal) {
						delete [] *iter;
            current->moves->erase(iter);
            return;
        }
    }
    assert(false); // If we got to here, then we did not find the point.
}

void
SupportTree::clear()
{
    delete root;
    root = new SupportNode();
}

const int32_t* SupportTree::reducable(
                    const int32_t* move,
                    const int32_t* skip_move) const
{
    //Statistics::incr_num_reducable_checks();
    return reducable(move, skip_move, root);
}

const int32_t* SupportTree::reducable_negative(
                    const int32_t* move,
                    const int32_t* skip_move) const
{
    //Statistics::incr_num_reducable_checks();
    return reducable_negative(move, skip_move, root);
}

/*
void
SupportTree::reducable(
                    const Binomial& b,
                    BinomialList& reducers) const
{
    reducable(b, reducers, root);
}
*/

const int32_t* SupportTree::reducable(
                    const int32_t* move,
                    const int32_t* skip_move,
                    const SupportNode* node) const
{
    assert(node != 0);

    for (int i = 0; i < (int) node->nodes.size(); ++i) {
        if (move[node->nodes[i].first] > 0) {
            const int32_t* reducer = reducable(move, skip_move, node->nodes[i].second);
            if (reducer != 0) {
								return reducer;
						}
        }
    }

    if (node->moves) {
        // Filter& f = *node->filter;
        for (MoveList::iterator i = node->moves->begin();
                        i != node->moves->end(); ++i) {
            const int32_t* cand_move = *i;
            if (reduces(cand_move, move)) {
                if (!skip_move || !equals(cand_move, skip_move)) { 
										return cand_move; 
								}
            }
        }
    }

    return 0;
}

const int32_t* SupportTree::reducable_negative(
                    const int32_t* move,
                    const int32_t* skip_move,
                    const SupportNode* node) const
{
    assert(node != 0);

    for (int i = 0; i < (int) node->nodes.size(); ++i) {
        if (move[node->nodes[i].first] < 0) {
            const int32_t* reducer = reducable_negative(move, skip_move, node->nodes[i].second);
            if (reducer != 0) 
						{ return reducer; }
        }
    }

    if (node->moves)
    {
        // Filter& f = *node->filter;
        for (MoveList::iterator i = node->moves->begin();
                        i != node->moves->end(); ++i) {
            const int32_t* cand_move = *i;
            if (reduces_negative(cand_move, move)) {
                if (skip_move && !equals(cand_move, skip_move)) {
										return cand_move;
								}
            }
        }
    }

    return 0;
}

/*
void
SupportTree::reducable(
                    const Binomial& b,
                    BinomialList& reducers,
                    const SupportNode* node) const
{
    assert(node != 0);

    for (int i = 0; i < (int) node->nodes.size(); ++i) {
        if (b[node->nodes[i].first] > 0) {
            reducable(b, reducers, node->nodes[i].second);
        }
    }

    if (node->binomials) {
        Filter& f = *node->filter;
        for (BinomialList::iterator i = node->binomials->begin();
                        i != node->binomials->end(); ++i) {
            const Binomial& bi = *(*i);
            if (Binomial::reduces(bi, f, b)) {
                reducers.push_back(&bi);
            }
        }
    }
}
*/


inline bool SupportTree::reduces(const int32_t* move1, const int32_t* move2) const
{
    for (int i = 0; i < dim; i++)
        if (move1[i] > 0 && move1[i] > move2[i])
            return false;
    return true;
}

inline bool SupportTree::reduces_negative(const int32_t* move1, const int32_t* move2) const
{
    for (int i = 0; i < dim; i++)
        if (move1[i] > 0 && move1[i] > -move2[i])
            return false;
    return true;

}

inline bool SupportTree::equals(const int32_t* move1, const int32_t* move2) const
{
    for (int i = 0; i < dim; i++)
        if (move1[i] != move2[i])
            return false;
    return true;
}

inline void SupportTree::print_move(const int32_t* move) const
{
    for (int i = 0; i < dim; i++)
				std::cout << move[i];
}


void SupportTree::print() const
{
		std::vector<int> labels;
    print(root, labels);
}

void SupportTree::print(const SupportNode* node, std::vector<int> &labels) const
{
    assert(node != 0);
    if (node->moves) {
        std::cout << "Num moves = " << node->moves->size() << std::endl;
				/*
        for (int i = 0; i < (int) node->filter->size(); ++i) {
            *out << (*node->filter)[i] << " ";
        }
				*/
        std::cout << "\n";
				std::cout << "[" << std::endl;
        for (MoveList::iterator i = node->moves->begin();
              i != node->moves->end(); ++i) {
						std::cout << "[";
						for (int j = 0; j < dim; j++)
            		std::cout << (*i)[j] << " ";
						std::cout << "]\n";
        }
				std::cout << "]";
				std::cout << std::endl;
    }
		else {
        std::cout << "Num moves = 0" << std::endl;
		}

    for (int i = 0; i < (int) node->nodes.size(); ++i) {
				labels.push_back(node->nodes[i].first);
				for (auto label : labels)
						std::cout << label << " ";
				std::cout << std::endl;
        print(node->nodes[i].second, labels);
				labels.pop_back();
    }
}

