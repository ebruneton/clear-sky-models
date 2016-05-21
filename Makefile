# Copyright (c) 2015 Eric Bruneton
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

GPP = g++
GPP_FLAGS = -Wall -Wmain -pedantic -pedantic-errors -std=c++11

DEBUG_FLAGS = -g
RELEASE_FLAGS = -DNDEBUG -O3 -fexpensive-optimizations

HEADERS := $(shell find . -name "*.h")
SOURCES := $(shell find . -name "*.cc" -not -name "*test*" -not -name "main.cc")
TEST_SOURCES := $(shell find . -name "*test*.cc")
ALL_SOURCES := $(HEADERS) $(SOURCES) $(TEST_SOURCES) main.cc
LINT_SOURCES := $(filter-out ./atmosphere/model/hosek/ArHosek%,$(ALL_SOURCES))

DEBUG_OBJECTS := $(SOURCES:%.cc=output/Debug/%.o)
RELEASE_OBJECTS := $(SOURCES:%.cc=output/Release/%.o)
TEST_OBJECTS := $(TEST_SOURCES:%.cc=output/Debug/%.o)

ARCHIVE_URL := \
    http://www.graphics.cornell.edu/resources/clearsky/data/2013-05-27/RADIANCE
ARCHIVES := \
    2013-05-27___09.30.00.7z 2013-05-27___09.45.00.7z \
    2013-05-27___10.00.00.7z 2013-05-27___10.15.00.7z \
    2013-05-27___10.30.00.7z 2013-05-27___10.45.00.7z \
    2013-05-27___11.00.00.7z 2013-05-27___11.15.00.7z \
    2013-05-27___11.30.00.7z 2013-05-27___11.45.00.7z \
    2013-05-27___12.00.00.7z 2013-05-27___12.15.00.7z \
    2013-05-27___12.30.00.7z 2013-05-27___12.45.00.7z \
    2013-05-27___13.00.00.7z 2013-05-27___13.15.00.7z \
    2013-05-27___13.30.00.7z

INPUTS := $(ARCHIVES:%.7z=input/%)

.DEFAULT = all
.PHONY = lint all test clean

output/Debug/%.o: %.cc
	mkdir -p $(@D)
	$(GPP) $(GPP_FLAGS) $(DEBUG_FLAGS) -I. -c $< -o $@

output/Release/%.o: %.cc
	mkdir -p $(@D)
	$(GPP) $(GPP_FLAGS) $(RELEASE_FLAGS) -I. -c $< -o $@

output/Debug/clearskymodels: $(DEBUG_OBJECTS) output/Debug/main.o
	mkdir -p $(@D)
	$(GPP) -s -o $@ $^

output/Release/clearskymodels: $(RELEASE_OBJECTS) output/Release/main.o
	mkdir -p $(@D)
	$(GPP) -s -o $@ $^

output/Debug/clearskymodels_test: $(DEBUG_OBJECTS) $(TEST_OBJECTS)
	mkdir -p $(@D)
	echo $(TEST_OBJECTS)
	$(GPP) -s -o $@ $^

lint: $(LINT_SOURCES)
	./cpplint.py --root=$(PWD) $^

input/%:
	wget -O input/$*.7z $(ARCHIVE_URL)/$*.7z
	7zr x -oinput input/$*.7z
	rm input/$*.7z

all: output/Debug/clearskymodels output/Release/clearskymodels $(INPUTS)
	mkdir -p output/cache/bruneton
	mkdir -p output/cache/haber
	mkdir -p output/cache/input
	mkdir -p output/cache/nishita
	mkdir -p output/cache/polradtran
	mkdir -p output/figures
	mkdir -p output/libradtran
	output/Release/clearskymodels \
	    /usr/local/bin/uvspec /usr/local/share/libRadtran/data

test: output/Debug/clearskymodels_test
	output/Debug/clearskymodels_test

clean:
	rm -rf output/Debug
	rm -rf output/Release

