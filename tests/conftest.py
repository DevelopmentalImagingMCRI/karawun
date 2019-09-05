# Copyright 2019 Murdoch Children's Research Institute,
# Melbourne, Australia
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--runcreatebaseline",
        action="store_true",
        default=False,
        help="create sha file"
    )


def pytest_configure(config):
    config.addinivalue_line("markers",
                            ("createbaseline: mark test as creating "
                             "baseline data"))


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runcreatebaseline"):
        # --runcreatebaseline given in cli: do not skip slow tests
        return
    skip_createbaseline = pytest.mark.skip(
        reason="need --runcreatebaseline option to run")
    for item in items:
        if "createbaseline" in item.keywords:
            item.add_marker(skip_createbaseline)
