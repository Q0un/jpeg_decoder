#include <fft.h>

#include <fftw3.h>
#include <cmath>

class DctCalculator::Impl {
public:
    Impl(size_t width, std::vector<double>* input, std::vector<double>* output)
        : input_(input), width_(width), output_(output) {
        if (width * width != input->size() || width * width != output->size()) {
            throw std::invalid_argument("Wrong input size");
        }
        plan_ = fftw_plan_r2r_2d(width, width, input->data(), output->data(), FFTW_REDFT01,
                                 FFTW_REDFT01, FFTW_MEASURE);
    }

    void Inverse() {
        for (size_t i = 0; i < width_; ++i) {
            input_->operator[](i) *= std::sqrt(2);
            input_->operator[](i * width_) *= std::sqrt(2);
        }
        fftw_execute(plan_);
        double inv = 1. / (width_ * 2);
        for (size_t i = 0; i < output_->size(); ++i) {
            output_->operator[](i) *= inv;
        }
    }

    ~Impl() {
        fftw_destroy_plan(plan_);
    }

private:
    fftw_plan plan_;
    std::vector<double>* input_;
    std::vector<double>* output_;
    size_t width_;
};

DctCalculator::DctCalculator(size_t width, std::vector<double>* input,
                             std::vector<double>* output) {
    impl_ = std::make_unique<Impl>(width, input, output);
}

void DctCalculator::Inverse() {
    impl_->Inverse();
}

DctCalculator::~DctCalculator() = default;
