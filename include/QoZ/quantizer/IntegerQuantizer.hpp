#ifndef _SZ_INTEGER_QUANTIZER_HPP
#define _SZ_INTEGER_QUANTIZER_HPP

#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include "QoZ/def.hpp"
#include "QoZ/quantizer/Quantizer.hpp"

namespace QoZ {

    template<class T>
    class LinearQuantizer : public concepts::QuantizerInterface<T> {
    public:
        LinearQuantizer() : error_bound(1), error_bound_reciprocal(1), radius(32768) {}

        LinearQuantizer(double eb, int r = 32768) : error_bound(eb),
                                                    error_bound_reciprocal(1.0 / eb),
                                                    radius(r) {
            assert(eb != 0);
        }

        int get_radius() const { return radius; }

        double get_eb() const { return error_bound; }

        void set_eb(double eb) {
            error_bound = eb;
            error_bound_reciprocal = 1.0 / eb;
        }

        void setTrimToZero(int ttz=1){
            trimToZero=ttz;
        }

        // quantize the data with a prediction value, and returns the quantization index
        int quantize(T data, T pred) {
            if(this->trimToZero and fabs(data)<=this->error_bound){
                
                return this->trimToZero-1;
            }
            T diff = data - pred;
            int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
            if (quant_index < this->radius * 2) {
                quant_index >>= 1;
                int half_index = quant_index;
                quant_index <<= 1;
                int quant_index_shifted;
                if (diff < 0) {
                    quant_index = -quant_index;
                    quant_index_shifted = this->radius - half_index;
                } else {
                    quant_index_shifted = this->radius + half_index;
                }
                T decompressed_data = pred + quant_index * this->error_bound;
                if (fabs(decompressed_data - data) > this->error_bound) {
                    return 0;
                } else {
                    return quant_index_shifted;
                }
            } else {
                return 0;
            }
        }

        // quantize the data with a prediction value, and returns the quantization index and the decompressed data
        // int quantize(T data, T pred, T& dec_data);
        int quantize_and_overwrite(T &data, T pred,bool save_unpred=true) {

            if(this->trimToZero>0 and fabs(data)<=this->error_bound){
                data=0;
                if(this->trimToZero==1 and save_unpred)
                    unpred.push_back(0);
                //std::cout<<data<<" "<<this->trimToZero-1<<std::endl;

                return this->trimToZero-1;
            }

            T diff = data - pred;
            int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
            if (quant_index < this->radius * 2) {
                quant_index >>= 1;
                int half_index = quant_index;
                quant_index <<= 1;
                int quant_index_shifted;
                if (diff < 0) {
                    quant_index = -quant_index;
                    quant_index_shifted = this->radius - half_index;
                } else {
                    quant_index_shifted = this->radius + half_index;
                }
                T decompressed_data = pred + quant_index * this->error_bound;
                if (fabs(decompressed_data - data) > this->error_bound or (this->trimToZero==2 and quant_index_shifted==1)) {
                    //std::cout<<data<<std::endl;
                    if(save_unpred)
                        unpred.push_back(data);
                    return 0;
                } else {
                    data = decompressed_data;
                    return quant_index_shifted;
                }
            } else {
                if(save_unpred)
                    unpred.push_back(data);
                return 0;
            }
        }

        int quantize_and_overwrite(T ori, T pred, T &dest,bool save_unpred=true) {

            if(this->trimToZero>0 and fabs(ori)<=this->error_bound){
                if(this->trimToZero==1 and save_unpred)
                    unpred.push_back(0);
                //std::cout<<ori<<" "<<this->trimToZero-1<<std::endl;
                dest = 0;
                return this->trimToZero-1;
            }
            T diff = ori - pred;
            int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
            if (quant_index < this->radius * 2) {
                quant_index >>= 1;
                int half_index = quant_index;
                quant_index <<= 1;
                int quant_index_shifted;
                if (diff < 0) {
                    quant_index = -quant_index;
                    quant_index_shifted = this->radius - half_index;
                } else {
                    quant_index_shifted = this->radius + half_index;
                }
                T decompressed_data = pred + quant_index * this->error_bound;
                if (fabs(decompressed_data - ori) > this->error_bound or (this->trimToZero==2 and quant_index_shifted==1)) {
                    if(save_unpred)
                        unpred.push_back(ori);
                    dest = ori;
                    return 0;
                } else {
                    dest = decompressed_data;
                    return quant_index_shifted;
                }
            } else {
                if(save_unpred)
                    unpred.push_back(ori);
                dest = ori;
                return 0;
            }
        }

        void insert_unpred(T ori){
            unpred.push_back(ori);
        }
        
        void print_unpred(){
            for(auto x:unpred)
                std::cout<<x<<std::endl;
        }

        // recover the data using the quantization index
        T recover(T pred, int quant_index) {

            if(this->trimToZero==2 and quant_index==1){
                return 0;
            }
            if (quant_index) {
                return recover_pred(pred, quant_index);
            } else {
                return recover_unpred();
            }
        }


        T recover_pred(T pred, int quant_index) {
            return pred + 2 * (quant_index - this->radius) * this->error_bound;
        }

        T recover_unpred() {
            return unpred[index++];
        }

        size_t size_est() {
            return unpred.size() * sizeof(T);
        }

        void save(unsigned char *&c) const {
            // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
            c[0] = 0b00000010;
            c += 1;
            // std::cout << "saving eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            *reinterpret_cast<double *>(c) = this->error_bound;
            
            c += sizeof(double);
            *reinterpret_cast<int *>(c) = this->radius;
           
            c += sizeof(int);
            *reinterpret_cast<size_t *>(c) = unpred.size();
           
            c += sizeof(size_t);
            memcpy(c, unpred.data(), unpred.size() * sizeof(T));
            c += unpred.size() * sizeof(T);
        };

        void load(const unsigned char *&c, size_t &remaining_length) {
            
            assert(remaining_length > (sizeof(uint8_t) + sizeof(T) + sizeof(int)));
            c += sizeof(uint8_t);
            remaining_length -= sizeof(uint8_t);
            this->error_bound = *reinterpret_cast<const double *>(c);
            //std::cout<<this->error_bound<<std::endl;
           
            this->error_bound_reciprocal = 1.0 / this->error_bound;
            c += sizeof(double);
            this->radius = *reinterpret_cast<const int *>(c);
            //std::cout<<this->radius<<std::endl;
            
            c += sizeof(int);
            size_t unpred_size = *reinterpret_cast<const size_t *>(c);
            //std::cout<<unpred_size<<std::endl;
            
            c += sizeof(size_t);
            this->unpred = std::vector<T>(reinterpret_cast<const T *>(c), reinterpret_cast<const T *>(c) + unpred_size);
            c += unpred_size * sizeof(T);
            // std::cout << "loading: eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            // reset index
            index = 0;
        }

        void print() {
            printf("[IntegerQuantizer] error_bound = %.8G, radius = %d, unpred = %lu\n", error_bound, radius, unpred.size());
        }

        void clear() {
            unpred.clear();
            index = 0;
        }


        virtual void postcompress_data() {
        }

        virtual void postdecompress_data() {
        }

        virtual void precompress_data() {};

        virtual void predecompress_data() {};


    private:
        std::vector<T> unpred;
        size_t index = 0; // used in decompression only

        double error_bound;
        double error_bound_reciprocal;
        int radius; // quantization interval radius
        int trimToZero=0;
    };

}
#endif
