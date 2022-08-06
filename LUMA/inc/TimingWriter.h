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

  TimingWriter(const std::string outdir, double dt)
    : outfile(outdir+"/timing.txt"),
      dt(dt),
      t_start(std::chrono::steady_clock::now()),
      t(0),
      runtime_s(0),
      speed(0),
      average_speed(0) {
    std::ofstream   os(outfile, std::ofstream::app);
    os << "timestep\tphysical_time\truntime[s]\taverage_speed_h[physical_time/h]" << std::endl;
    writeTimingData();
  }

  void recordTimingData(int t) {
    double previous_runtime_s = runtime_s;
    int previous_t = this->t;
    auto t_now = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = t_now - t_start;

    this->t = t;
    runtime_s = duration.count();
    speed = (runtime_s - previous_runtime_s) / ((t - previous_t) * dt);
    average_speed = (t * dt) / runtime_s;
  }

  void writeTimingData() {
    std::ofstream   os(outfile, std::ofstream::app);
    os << t << "\t" << t * dt << "\t" << runtime_s << "\t" << average_speed*3600 << std::endl;
  }

private:
  std::string outfile;
  double dt;
  std::chrono::time_point<std::chrono::steady_clock> t_start;

  int t;
  double runtime_s;
  double speed;
  double average_speed;
};

#endif
