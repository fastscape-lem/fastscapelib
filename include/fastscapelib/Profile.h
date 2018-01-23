/**
Copyright (c) 2017, INPG
Main authors : Guillaume Cordonnier, Antoine Begault
All rights reserved.

Redistribution and use in source and binary forms of this module (see bellow for a complete list of what is not included in this license), with or without modification, are permitted only for non-commercial uses, provided that the following conditions are met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
* Neither the name of the INPG nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE REGENTS AND CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef PROFILE_H
#define PROFILE_H

#include <chrono>
#include <string>
#include <iostream>

// minimal profiler. Each "count" (default = 10) call,
// it gives the min, max, and average time of execution of the current scope

// it needs a unique type identifier as "type_name".

// exemples:

/*

  void f()
  {
    PROFILE(profile_f, "f");

    // do something
  }

  void g()
  {
    // ..

    {
      PROFILE_COUNT(scope_in_g, "Cool stuff", 100);
      // cool stuff
    }

    //..
  }

*/



#define PROFILE(type_name, name) class type_name; Profile<type_name> profile_scope_var_ ## type_name ## _sjflfo (name)
#define PROFILE_COUNT(type_name, name, count) class type_name; Profile<type_name, count> profile_scope_var_ ## type_name ## count ## _sjflfo (name)

namespace ns_profile_types {


class StdTimer
{
public:
    StdTimer(){}

    void start()
    {
        chrono_start_  = std::chrono::high_resolution_clock::now();
    }
    double get_ms()
    {
        chrono_stop_ = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = chrono_stop_-chrono_start_;

        return elapsed_seconds.count() *1000.0;
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> chrono_start_, chrono_stop_;
};

}


template<typename id, class Timer, int count>
class Profile_t
{
public:
    Profile_t(std::string name = "") :
        name_ {name}
    {
        timer.start();
    }

    ~Profile_t()
    {
        double mseconds = timer.get_ms();

                min_ = cur_count ? std::min(min_, mseconds) : mseconds;
        max_ = cur_count ? std::max(max_, mseconds) : mseconds;
        sum_ = cur_count ? sum_ + mseconds : mseconds;

        if(++cur_count == count)
        {
            /*std::cerr << "Profile " << name_
                      << ": min = " << min_
                      << "ms, max = " << max_
                      << "ms, avg = " << sum_ / double(count)
                      << "ms" << std::endl;*/
            cur_count = 0;
        }

    }

private:
    std::string name_;

    Timer timer;

    static int cur_count;
    static double min_, max_ , sum_;

};

template<typename id, class Timer, int count>
int Profile_t<id, Timer, count>::cur_count = 0;

template<typename id, class Timer, int count>
double Profile_t<id, Timer, count>::min_ = 0.0;

template<typename id, class Timer, int count>
double Profile_t<id, Timer, count>::max_ = 0.0;

template<typename id, class Timer, int count>
double Profile_t<id, Timer, count>::sum_ = 0.0;


template <typename id, int count = 3>
using Profile = Profile_t<id, ns_profile_types::StdTimer, count>;

#endif // PROFILE_H
