#ifndef UNION_FIND_TREE_HPP__
#define UNION_FIND_TREE_HPP__

#include <vector>

namespace openjij {
    namespace utility {
        struct UnionFindTree {
            using Parent = std::vector<size_t>;
            using Rank = std::vector<size_t>;

            UnionFindTree(const size_t n);

            void unite(const size_t x, const size_t y);
            int get_root(const size_t x) const;
            bool is_same_root(const size_t x, const size_t y) const;

        private:
            Parent _parent;
            Rank _rank;
        };
    } // namespace utility
} // namespace openjij

#endif /* end of include guard */
