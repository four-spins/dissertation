#include <types.hpp>
#include <vector>
#include <cmath>
/// statistics provides some convenience functions 
/// -- jackknife algorithm
namespace statistics
{
    /// jackknife is a class that provides a jackknife algorithm
    /// for correct error estimates on data of type double\n
    /// \n
    /// usage:\n
    /// 
    /// double f(std::vector<doube>  &) { calculate a double };\n
    /// statistics::jackknife foo(data, block size, f);\n
    /// double jk_mean = foo.get_mean();\n
    /// double jk_std = foo.get_standard_deviation();\n
    ///
    /// @warning Be careful, how the function is implemented!
    /// if you want, for example calculate an estimator for x**2
    /// then provide:\n
    ///    double f(std::vector<double> &x) { return mean(x)**2 }\n
    /// and NOT:\n
    ///    double f(std::vector<double> &x) { return mean(x**2) }\n
    /// The latter function would simply calculate the bad estimator
    /// that we try to avoid using jackknife algorithm !
    ///
    class jackknife
    {
        public: 
            jackknife
            (
                std::vector<double> &data_,
                uint32 blocksize_,
                double (*f_)(std::vector<double> &)
            );
            double get_mean();
            double get_standard_deviation();

        private:
            // counts how many blocks were already calculated
            uint32 nr_blocks;
            // number of elements in one block
            uint32 blocksize;
            // total number of datapoints
            uint32 nr_datapoints;
            // complete data
            std::vector<double> &data;
            // Blocks
            std::vector<double> block_values;
            double jk_mean, jk_std;

            // (arbitrary) function to be executed by jackknife algorithm
            double (*f)(std::vector<double> &);
            void calculate_jackknife_blocks();
    };

    // constructor
    // for cases where we stream data into the algorithm
    jackknife::jackknife
    (
        std::vector<double> &data_,
        uint32 blocksize_,
        double (*f_)(std::vector<double> &)
    )
    : data(data_), blocksize(blocksize_), f(f_)
    {
        nr_blocks = data.size()/blocksize;
        if (nr_blocks*blocksize != data.size())
            std::cerr << "Warning: Jackknife-algorithm skips some data points"
                      << " as blocksize and data.size() don't match\n";
        
        // allocate some memory for performance
        // we don't know the actual amount of 
        block_values.resize(nr_blocks);
        calculate_jackknife_blocks();
    }


    void jackknife::calculate_jackknife_blocks()
    {

        std::cout << "data: \n";
        for (int i = 0; i< data.size(); i++)
        {
            std::cout << data[i] << "\t";
        }

        std::vector<double> temp_block(blocksize, 0);
        double f_all = f(data);
        double sum = 0.0;
        for(uint32 i = 0; i < nr_blocks; i++)
        {
            for (uint32 j=0; j < blocksize; j++)
            {
                temp_block.at(j) = data.at(i*blocksize + j);
            }
            block_values.at(i) = (data.size()*f_all - blocksize*f(temp_block))
                /static_cast<double>(data.size() - blocksize);
            sum += block_values.at(i);
        }
        
        jk_mean = sum/static_cast<double>(block_values.size());

        sum = 0.0;
        for (uint32 i = 0; i < nr_blocks; i++)
            sum += pow((block_values.at(i) - jk_mean),2);
        
        std::cout << "Nr of blocks: " << nr_blocks << "\n";
        double foo = static_cast<double>(nr_blocks);
        jk_std = sqrt((foo-1.0)/foo * sum);

    }

    double jackknife::get_mean()
    {
        return jk_mean;
    }

    double jackknife::get_standard_deviation()
    {
        return jk_std;    
    }

    class discrete_histogram
    {
        public:
            discrete_histogram(uint32 nr_bins_, uint32 value_min_);
            void add(uint32 value);
            uint64 get(uint32 value);

        private:
            std::vector<uint64> histogram;
            uint32 nr_bins;
            uint32 value_min;
    };

    discrete_histogram::discrete_histogram(uint32 nr_bins_, uint32 value_min_)
    : nr_bins(nr_bins_), value_min(value_min_)
    {
        histogram.resize(nr_bins);
        for (uint32 i=0; i<nr_bins; i++)
        {
            histogram.at(i) = 0;
        }
    }
    
    void discrete_histogram::add(uint32 value)
    {
        uint32 index = value - value_min;
        histogram.at(index) ++;
    }

    uint64 discrete_histogram::get(uint32 value)
    {
        return histogram.at(value - value_min);    
    }
}
