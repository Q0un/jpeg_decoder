#include "huffman.h"

#include <memory>

namespace {

struct Node {
    Node() = default;
    virtual ~Node() = default;

    std::shared_ptr<Node> left;
    std::shared_ptr<Node> right;
};

struct Leaf : public Node {
    Leaf(uint8_t value) : value(value) {
    }

    uint8_t value;
};

}  // anonymous namespace

class HuffmanTree::Impl {
public:
    void Build(std::vector<uint8_t> code_lengths, const std::vector<uint8_t>& values) {
        if (code_lengths.size() > 16) {
            throw std::invalid_argument("Maximum size of code length is 16");
        }
        root_ = std::make_unique<Node>();
        cur_ = root_;
        size_t last = 0;
        RecursiveBuild(root_, 0, last, code_lengths, values);
        if (last < values.size()) {
            root_ = cur_ = nullptr;
            throw std::invalid_argument("Can't build huffman tree with this args");
        }
        for (const auto& i : code_lengths) {
            if (i != 0) {
                root_ = cur_ = nullptr;
                throw std::invalid_argument("Can't build huffman tree with this args");
            }
        }
    }

    bool Move(bool bit, int& value) {
        if (!cur_) {
            throw std::invalid_argument("Huffman tree is empty");
        }
        if (bit) {
            cur_ = cur_->right;
        } else {
            cur_ = cur_->left;
        }
        if (const auto p = std::dynamic_pointer_cast<Leaf>(cur_)) {
            value = p->value;
            cur_ = root_;
            return true;
        }
        return false;
    }

private:
    void RecursiveBuild(std::shared_ptr<Node> node, size_t depth, size_t& last,
                        std::vector<uint8_t>& code_lengths, const std::vector<uint8_t>& values) {
        if (last == values.size()) {
            return;
        }
        if (code_lengths[depth]) {
            node->left = std::make_unique<Leaf>(values[last++]);
            --code_lengths[depth];
        } else {
            node->left = std::make_unique<Node>();
            RecursiveBuild(node->left, depth + 1, last, code_lengths, values);
        }
        if (last == values.size()) {
            return;
        }
        if (code_lengths[depth]) {
            node->right = std::make_unique<Leaf>(values[last++]);
            --code_lengths[depth];
        } else {
            node->right = std::make_unique<Node>();
            RecursiveBuild(node->right, depth + 1, last, code_lengths, values);
        }
    }

    std::shared_ptr<Node> root_;
    std::shared_ptr<Node> cur_;
};

HuffmanTree::HuffmanTree() {
    impl_ = std::make_unique<HuffmanTree::Impl>();
}

void HuffmanTree::Build(const std::vector<uint8_t>& code_lengths,
                        const std::vector<uint8_t>& values) {
    impl_->Build(code_lengths, values);
}

bool HuffmanTree::Move(bool bit, int& value) {
    return impl_->Move(bit, value);
}

HuffmanTree::HuffmanTree(HuffmanTree&&) = default;

HuffmanTree& HuffmanTree::operator=(HuffmanTree&&) = default;

HuffmanTree::~HuffmanTree() = default;
