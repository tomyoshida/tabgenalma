# tabgenalma
Generate a latex table of an observation summary from ALMA measurement sets.

# Usage
```python
from tabgenalma import generate_table

vis_list = [ 'band3.ms', 'band6.ms' ]
bands = [ 3, 6 ]

generate_table( vis_list, bands )
```
