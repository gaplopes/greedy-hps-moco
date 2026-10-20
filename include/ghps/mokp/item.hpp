#ifndef MOKP_ITEM_HPP
#define MOKP_ITEM_HPP

#include <cstdint>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace moco::ghps {


struct Item {
 public:
  int32_t idx;                  // Original input row index
  int64_t weight;               // Weight of the item
  std::vector<int64_t> values;  // Values of the item for each objective

  Item(int idx, int64_t weight, std::vector<int64_t> values)
      : idx(idx), weight(weight), values(std::move(values)) {}

  std::string to_string() const {
    std::stringstream ss;
    ss << "Item(" << idx << ", " << weight << ", [";
    for (int i = 0; i < (int)values.size(); i++) {
      ss << values[i];
      if (i < (int)values.size() - 1) {
        ss << ", ";
      }
    }
    ss << "])\n";
    return ss.str();
  }
};

} // namespace moco::ghps

#endif // MOKP_ITEM_HPP
