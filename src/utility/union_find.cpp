#include "union_find.h"

namespace openjij {
    namespace utility {
        UnionFindTree::UnionFindTree(const size_t n)
            : _parent(n), _rank(n) {
            for (size_t i = 0; i < n; ++i) {
                _parent[i] = i;
            }
        }

        void UnionFindTree::unite(const size_t x, const size_t y) {
            auto root_x = get_root(x);
            auto root_y = get_root(y);

            if (root_x == root_y) {
                return;
            }

            if (_rank[root_x] < _rank[root_y]) {
                _parent[root_x] = root_y;
            } else {
                _parent[root_y] = root_x;

                if (_rank[root_x] == _rank[root_y]) {
                    ++_rank[root_x];
                }
            }
        }

        int UnionFindTree::get_root(const size_t x) const {
            if (_parent[x] == x) {
                return x;
            } else{
                return get_root(_parent[x]);
            }
        }

        bool UnionFindTree::is_same_root(const size_t x, const size_t y) const {
            auto root_x = get_root(x);
            auto root_y = get_root(y);
            return root_x == root_y;
        }
    } // namespace utility
} // namespace openjij
