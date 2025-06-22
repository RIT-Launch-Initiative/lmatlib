# Launch MATLAB libraries
MATLAB utilities for data import/export from simulations, clean plotting, and simulation.

## Layout
```
data/     - data import/export
plot/     - plotting
sim/      - simulation
samples/  - usage examples
```

## How to download

Click `New->Project->From Git` and link to the repository (the last part should be `lmatlib`, with no `/tree/main`). This pulls all required files and lets you pull new changes as they are published.

## How to use

If you only want to run `samples`, just open the project file `Lmatlib.prj`.

If you want to use this as part of a different project, [Create a MATLAB Project](https://www.mathworks.com/help/matlab/matlab_prog/create-projects.html) under a common folder with `lmatlib`, then create a [*relative* reference](https://www.mathworks.com/help/simulink/ug/add-or-remove-a-reference-to-another-project.html) `lmatlib`. If projects are always created in this way, there only needs to be one `lmatlib` on your system that every other project references. The resulting folder structure should look like this:
```
<Common_folder_name>/
   <new-project>/ 
        <new-project-file>.prj
        [stuff]
   lmatlib/
        Lmatlib.prj
        [stuff]
```


## Contribution

Help line format---replace text inside and including angle brackets `<>`
```
%% <One-liner description>
% [<output1>, <output2>] = <function name>(<input1>, <input2>, <name1 = value1>, <name2 = value2>)
% Inputs
%   <input1>     (<data type(s) or allowed values>)  description
%   <input2>     (<data type(s) or allowed values>)  description
% <Name-value inputs, if present -- stuff like DisplayName = "blarg">
%   <name1>      (<data type(s) or allowed values>)  description
%   <name2>      (<data type(s) or allowed values>)  description
% Outputs
%   <output1>    (<data type(s) or allowed values>)  description
%   <output2>    (<data type(s) or allowed values>)  description
%
% EXAMPLES
% % no name-value arguments behaves one way
% [a, b] = func(in1, in2)
% % 
% 
```
