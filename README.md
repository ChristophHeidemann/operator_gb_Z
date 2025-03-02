# oprator_gb_Z

## Description
This package is a modification of https://github.com/ClemensHofstadler/operator_gb that allows for Gröbner basis computations with integer coefficients.

## License
Distributed under the terms of the GNU General Public License (GPL, see the LICENSE file), either version 2 or (at your option) any later version
- http://www.gnu.org/licenses/

## Technologies Used
- **Programming Languages**: SageMath
- **Libraries/Frameworks**: pyahocorasick

## Requirements:
- Sagemath 9.1 or later

## Setup Instructions:
Clone the repository:
   ```bash
   git clone https://github.com/ChristophHeidemann/operator_gb_Z.git
````

## Example:
```python
from operator_gb_Z import *
F.<x,y> = FreeAlgebra(ZZ,2)               
A = MyFreeAlgebra(ZZ,F.gens())
I = NCIdeal([2*x*y*x + x, 3*y*x- y])
I.groebner_basis(4)
```
