/**
 Copyright (c) 2015 Eric Bruneton
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
// Based on "Computing the Solar Vector" by Manuel Blanco-Muriel,
// Diego C. Alarcon-Padilla, Teodoro Lopez-Moratalla, and Martin Lara-Coira,
// in "Solar energy", vol 27, number 5, 2001 by Pergamon Press.
// Original code available at http://www.psa.es/sdg/sunpos.htm.

#include "atmosphere/sun_direction.h"

#include <cmath>
#include <cstdint>

// Declaration of some constants.
#define PI 3.14159265358979323846
#define TWO_PI (2 * PI)
#define RAD (PI / 180)
#define EARTH_MEAN_RADIUS 6371.01    // In km
#define ASTRONOMICAL_UNIT 149597890  // In km

void GetSunDirection(const Time& udt_time, const Location& udt_location,
    SunCoordinates *udt_sun_coordinates) {
  // Main variables.
  double elapsed_julian_days;
  double decimal_hours;
  double ecliptic_longitude;
  double ecliptic_obliquity;
  double right_ascension;
  double declination;

  // Auxiliary variables.
  double dy;
  double dx;

  // Calculate difference in days between the current Julian Day
  // and JD 2451545.0, which is noon 1 January 2000 Universal Time.
  {
    double julian_date;
    int64_t aux1;
    int64_t aux2;
    // Calculate time of the day in UT decimal hours.
    decimal_hours =
        udt_time.hours + (udt_time.minutes + udt_time.seconds / 60.0) / 60.0;
    // Calculate current Julian Day.
    aux1 = (udt_time.month - 14) / 12;
    aux2 = (1461 * (udt_time.year + 4800 + aux1)) / 4 +
        (367 * (udt_time.month - 2 - 12 * aux1)) / 12 -
        (3 * ((udt_time.year + 4900  + aux1) / 100)) / 4 +
        udt_time.day - 32075;
    julian_date = static_cast<double>(aux2) - 0.5 + decimal_hours / 24.0;
    // Calculate difference between current Julian Day and JD 2451545.0
    elapsed_julian_days = julian_date - 2451545.0;
  }

  // Calculate ecliptic coordinates (ecliptic longitude and obliquity of the
  // ecliptic in radians but without limiting the angle to be less than 2*Pi
  // (i.e., the result may be greater than 2*Pi).
  {
    double mean_longitude;  // Radians
    double mean_anomaly;
    double omega;
    omega = 2.1429 - 0.0010394594 * elapsed_julian_days;
    mean_longitude = 4.8950630 + 0.017202791698 * elapsed_julian_days;
    mean_anomaly = 6.2400600 + 0.0172019699 * elapsed_julian_days;
    ecliptic_longitude = mean_longitude + 0.03341607 * sin(mean_anomaly)  +
        0.00034894 * sin(2 * mean_anomaly) - 0.0001134 -
        0.0000203 * sin(omega);
    ecliptic_obliquity = 0.4090928 - 6.2140e-9 * elapsed_julian_days +
        0.0000396 * cos(omega);
  }

  // Calculate celestial coordinates (right ascension and declination) in
  // radians but without limiting the angle to be less than 2*Pi (i.e., the
  // result may be greater than 2*Pi).
  {
    double sin_ecliptic_longitude;
    sin_ecliptic_longitude =  sin(ecliptic_longitude);
    dy = cos(ecliptic_obliquity) * sin_ecliptic_longitude;
    dx = cos(ecliptic_longitude);
    right_ascension = atan2(dy, dx);
    if (right_ascension < 0.0) {
      right_ascension = right_ascension + TWO_PI;
    }
    declination = asin(sin(ecliptic_obliquity) * sin_ecliptic_longitude);
  }

  // Calculate local coordinates (azimuth and zenith angle) in degrees.
  {
    double greenwich_mean_sidereal_time;
    double local_mean_sidereal_time;
    double latitude_in_radians;
    double hour_angle;
    double cos_latitude;
    double sin_latitude;
    double cos_hour_angle;
    double parallax;
    greenwich_mean_sidereal_time =
        6.6974243242 + 0.0657098283 * elapsed_julian_days + decimal_hours;
    local_mean_sidereal_time =
        (greenwich_mean_sidereal_time * 15 + udt_location.longitude) * RAD;
    hour_angle = local_mean_sidereal_time - right_ascension;
    latitude_in_radians = udt_location.latitude * RAD;
    cos_latitude = cos(latitude_in_radians);
    sin_latitude = sin(latitude_in_radians);
    cos_hour_angle = cos(hour_angle);
    udt_sun_coordinates->zenith = acos(cos_latitude * cos_hour_angle
      * cos(declination) + sin(declination) * sin_latitude);
    dy = -sin(hour_angle);
    dx = tan(declination) * cos_latitude - sin_latitude  *cos_hour_angle;
    udt_sun_coordinates->azimuth = atan2(dy, dx);
    if (udt_sun_coordinates->azimuth < 0.0) {
      udt_sun_coordinates->azimuth = udt_sun_coordinates->azimuth + TWO_PI;
    }
    udt_sun_coordinates->azimuth = udt_sun_coordinates->azimuth / RAD;
    // Parallax Correction.
    parallax = (EARTH_MEAN_RADIUS / ASTRONOMICAL_UNIT) *
        sin(udt_sun_coordinates->zenith);
    udt_sun_coordinates->zenith =
        (udt_sun_coordinates->zenith  + parallax) / RAD;
  }
}
