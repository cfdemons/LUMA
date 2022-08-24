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
#include <limits>
#include <vector>

struct TimingPoint
{
	int t;
	double runtime;
};

class TimingHistory
{
public:

	TimingHistory(int len) : len(len) {
		points.resize(len);
	}

	void addPoint(const TimingPoint &point) {
		if (used == 0) {
			iNewest = 0;
			iOldest = 0;
			points.at(iNewest) = point;
			used++;
		} else if (used < len) {
			iNewest++;
			points.at(iNewest) = point;
			used++;
		} else if (used == len) {
			iNewest = (iNewest + 1) % len;
			iOldest = (iOldest + 1) % len;
			points.at(iNewest) = point;
		} else {
			throw std::range_error("TimingHistory::addPoint: used out of range");
		}
	}

	TimingPoint &getNewest() {
		assert(used > 0);
		return points.at(iNewest);
	}

	TimingPoint &getOldest() {
		assert(used > 0);
		return points.at(iOldest);
	}

private:
	int iNewest = -1;
	int iOldest = -1;
	int len;
	int used = 0;
	std::vector<TimingPoint> points;
};



class TimingWriter
{
public:

	TimingWriter(const std::string &outdir, double dt, int total_timesteps, int averaging_window,
				 bool write)
		: outfile(outdir+"/timing.txt"),
		  dt(dt),
		  total_timesteps(total_timesteps),
		  t_start(std::chrono::steady_clock::now()),
		  t(0),
		  runtime_s(0),
		  speed(0),
		  average_speed(0),
		  moving_average_speed(0),
		  runtime_remaining_h(std::numeric_limits<double>::quiet_NaN()),
		  history(averaging_window) {
		if (write) {
			std::ofstream   os(outfile, std::ofstream::app);
			os << "timestep\tphysical_time\truntime[s]\tspeed[physical_time/h]\taverage_speed[physical_time/h]\tmoving_average_speed[physical_time/h]\truntime_remaining[h]" << std::endl;
			writeTimingData();
		}
		history.addPoint(TimingPoint{0, 0});
	}

	void recordTimingData(int t) {
		double previous_runtime_s = runtime_s;
		int previous_t = this->t;
		auto t_now = std::chrono::steady_clock::now();
		std::chrono::duration<double> duration = t_now - t_start;

		this->t = t;
		runtime_s = duration.count();
		history.addPoint(TimingPoint{t, runtime_s});

		speed = ((t - previous_t) * dt) / (runtime_s - previous_runtime_s);
		average_speed = (t * dt) / runtime_s;
	
		TimingPoint hNewest = history.getNewest();
		TimingPoint hOldest = history.getOldest();

		moving_average_speed = dt * (hNewest.t - hOldest.t) /
			(hNewest.runtime - hOldest.runtime);

		runtime_remaining_h = (total_timesteps - t) * dt / moving_average_speed / 3600;
	}

	void writeTimingData() {
		std::ofstream   os(outfile, std::ofstream::app);
		os << 
			t << "\t" << 
			t * dt << "\t" << 
			runtime_s << "\t" << 
			speed*3600 << "\t" <<
			average_speed*3600 << "\t" <<
			moving_average_speed*3600 << "\t" <<
			runtime_remaining_h <<
			std::endl;
	}

private:
	std::string outfile;
	double dt;
	int total_timesteps;
	std::chrono::time_point<std::chrono::steady_clock> t_start;

	int t;
	double runtime_s;
	double speed;
	double average_speed;
	double moving_average_speed;
	double runtime_remaining_h;

	TimingHistory history;
};

#endif
