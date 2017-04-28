/******************************************************************************
 *
 * @file: timer.hpp
 *
 * @date: 09/22/2012 01:18:09 PM (CEST)
 *
 * @author: Marco Müller <muelma@gmail.com>
 *
 ******************************************************************************/
# ifndef TIMER_HPP 
# define TIMER_HPP 
# include <chrono>

// use it like
//
//    timer estim_time;
//    // do something
//    unsigned seconds = estim_time.elapsed().count();
//
class timer
{   
    private:
        typedef std::chrono::high_resolution_clock clock;
        clock::time_point start;

    public:
        timer() : start(clock::now()){}
        clock::time_point restart() { start = clock::now(); return start;}
        std::chrono::milliseconds elapsed_ms()   
            { return std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - start);}
        std::chrono::seconds elapsed()
            { return std::chrono::duration_cast<std::chrono::seconds>(clock::now() - start);}
};
# endif // TIMER_HPP
