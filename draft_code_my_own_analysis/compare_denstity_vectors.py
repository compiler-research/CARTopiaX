# -----------------------------------------------------------------------------
# Copyright (C) 2025 Salvador de la Torre Gonzalez
# Co-author: Luciana Melina Luque
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -----------------------------------------------------------------------------

import numpy as np

with open('vector_densities.csv', 'r') as file:
    line = file.readline().strip()
    vector_theirs = np.array([float(val) for val in line.split(',')])

print(f"Vector length: {len(vector_theirs)}")

with open('vector_densities_mine.csv', 'r') as file:
    line = file.readline().strip()
    vector_mine = np.array([float(val) for val in line.split(',')])

print(f"Vector length: {len(vector_mine)}")


if np.array_equal(vector_theirs, vector_mine):
    print("The vectors are equal.")
else:
    print("The vectors are different at the following positions:")
    for i, (a, b) in enumerate(zip(vector_theirs, vector_mine)):
        if a != b:
            print(f"Index {i}: vector_theirs={a}, vector_mine={b}")