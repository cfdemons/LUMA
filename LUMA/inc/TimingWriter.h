/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
* Copyright 2022 The University of Manchester
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.*
*/

#ifndef TIMING_H
#define TIMING_H

#include <chrono>
#include <fstream>

class TimingWriter
{
public:

  TimingWriter(const std::string outdir)
    : t_start(std::chrono::steady_clock::now()),
      t(-1),
      runtime_s(-1),
      speed(-1),
      outfile(outdir+"/timing.txt") {
    std::ofstream   os(outfile, std::ofstream::app);
    os << "timestep\truntime[s]\tspeed_s[timesteps/s]\taverage_speed_s[timesteps/s]\taverage_speed_h[timesteps/h]" << std::endl;
    writeTimingData(0, 0.0, 0.0, 0.0, 0.0);
  }

  void recordTimingData(int t) {
    double previous_runtime_s = runtime_s;
    int previous_t = this->t;
    auto t_now = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = t_now - t_start;
    runtime_s = duration.count();
    this->t = t;
    speed = (runtime_s - previous_runtime_s) / (t - previous_t);
    average_speed = runtime_s / t;
  }

  void writeTimingData() {
    writeTimingData(t, runtime_s, speed, average_speed, average_speed * 3600);
  }

private:

  void writeTimingData(int t, double runtime_s, double speed, double average_speed) {
    std::ofstream   os(outfile, std::ofstream::app);
    os << t << "\t" << runtime_s << "\t" << speed << "\t" << average_speed << std::endl;
  }

  std::chrono::time_point<std::chrono::steady_clock> t_start;
  int t;
  double runtime_s;
  double speed;
  double average_speed;
  std::string outfile;
  std::vector<double> runtime_history;
};

#endif
