# Concentration vs. Time Graphing
Make Concentration vs. Time graphs for multi-molecular and multi-step reactions

Input the reactions as strings into a matrix, the rate constants into a matrix, and the initial concentration of each species into a dictionary

Will automatically detect when the reactions have reached steady-state and will adjust abscissa accordingly

Example: 

reactions = ['A+B=C', 'B+C=D']

ks = [1, 1.5]

C0 = {'A': 1, 'B': 1, 'C': 0, 'D': 0}

reactiongraphing(reactions, ks, C0)

![image](https://github.com/user-attachments/assets/6e22fbff-c626-4f45-aa86-936b54efce94)

