#include "../include/fitness.hpp"
#include <stdexcept>

using namespace EvoAlg;

size_t AbstractFitness::getSize() const {
    return size_;
}

int8_t AbstractFitness::getSign(size_t index) const {
    if (index >= size_) {
        throw std::out_of_range("out of range access");
    }
    return sign_[index];
}