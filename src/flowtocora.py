import sys


class LinHybridFlowStarToCORA:

    def __readFile(self, infile, options, name):
        """
        This method read a Flow* file and converts every part of it into CORA code.
        :param infile: The Flow* file which has to be converted
        :param options: Options relevant for CORA
        :param name: The name of the function and file in CORA
        :return: string containing the CORA code
        """
        count = 0
        init_line = 0
        modes_line = 0
        jumps_line = 0
        unsafe_line = 0

        with open(infile, 'r') as file:
            for line in file:
                count += 1
                if 'state var ' in line:
                    line = line.replace('state var ', '')
                    line = line.replace('\n', '')
                    line = line.replace(' ', '')
                    state_var = line.split(',')
                    options['stateVars'] = state_var
                if 'fixed steps' in line:
                    stepSize = line.split(' ')[4]
                    stepSize = stepSize.replace('\n', '')
                    options['stepSize'] = stepSize
                if 'time' in line:
                    tFinal = line.split(' ')[3]
                    tFinal = tFinal.replace('\n', '')
                    options['tFinal'] = tFinal
                if 'remainder estimation' in line:
                    reminder = line.split(' ')[4]
                    reminder = reminder.replace('\n', '')
                    options['reminder'] = reminder
                if 'max jumps' in line:
                    jumps = line.split(' ')[4]
                    jumps = jumps.replace('\n', '')
                    options['jumps'] = jumps
                if 'init' in line:
                    init_line = count
                if 'modes' in line:
                    modes_line = count
                if 'jumps' in line and 'max' not in line:
                    jumps_line = count
                if 'unsafe' in line:
                    unsafe_line = count

        file.close()

        res = ""
        res += self.__getInitialization(init_line, infile, name, options)
        res += self.__writeOptions(options)
        temp, loc_names, inv_names = self.__getModes(modes_line, infile, options)
        res += temp
        res += self.__getJumps(jumps_line, infile, options['stateVars'], loc_names)
        res += self.__defineLocations(loc_names, inv_names)
        res += self.__defineHybridAutomaton()
        res += self.__drawReachableSet()

        return res

    def __getInitialization(self, start, infile, name, options):
        res = "function res = " + name + "()\n\n"
        res += "%set options---------------------------------------------------------------\n"
        initial = []

        x0 = "options.x0 = ["
        r0 = "options.R0 = zonotope([options.x0, "

        with open(infile, 'r') as file:
            for i in range(start + 2):
                file.readline()
                if i == start:
                    initState = file.readline()
                    initState = initState.replace('\n', '')
                    initState = initState.replace(' ', '')
                    options['initState'] = initState

            # Now, extract the initalization information
            for line in file:
                if '}' in line:
                    break
                else:
                    l = line.split(' ')
                    initial.append(l[len(l) - 1])
        file.close()

        count = 0
        for item in initial:
            item = item.replace('[', '')
            item = item.replace(']', '')
            limits = item.split(',')
            lower = float(limits[0])
            upper = float(limits[1])
            middle = (upper - lower) / 2
            x0 += str(lower + middle)
            if count < len(initial) - 1:
                x0 += "; "
            count += 1

        x0 += "];\n"
        r0 += options['stepSize'] + " * eye(" + str(len(initial)) + ")]);\n"

        res += x0
        res += r0
        return res

    def __getModes(self, start, infile, options):
        """
        This method extracts the locations from the Flow* file.
        :param start: Line at which the modes specifications start
        :param infile: The input file
        :param options: Relevant options for CORA
        :return: string containing all the important information for locations
        """
        res = ""
        inv_area = False
        loc_names = []
        flows = []
        open_braces = 1
        invariants = []

        with open(infile, 'r') as file:
            for i in range(start + 1):
                file.readline()
            for line in file:
                if open_braces == 0:
                    break
                if '{' in line:
                    open_braces += 1
                if '}' in line:
                    open_braces -= 1
                if 'inv' in line:
                    inv_area = True
                    line = file.readline()
                if '{' in line and open_braces == 3 and inv_area == False:
                    # Save flow
                    fl = []
                    line = file.readline()
                    while '}' not in line:
                        f = line.replace('\n', '')
                        fl.append(f)
                        line = file.readline()
                    flows.append(fl)
                if '{' in line and open_braces == 3 and inv_area:
                    # Save invariant
                    inv = []
                    line = file.readline()
                    while '}' not in line:
                        it = line.replace('\n', '')
                        if 'in' in it:
                            # Invariant is given in interval representation
                            line_array = it.split('in')
                            var = line_array[0]
                            interval = line_array[1].replace('[','')
                            interval = interval.replace(']', '')
                            interval = interval.split(',')
                            inv.append(var + " >= " + interval[0])
                            inv.append(var + " <= " + interval[1])
                        else:
                            inv.append(it)
                        line = file.readline()
                    invariants.append(inv)
                    inv_area = False
                    open_braces -= 1
                if open_braces == 1 and '}' not in line:
                    line = line.replace(' ', '')
                    line = line.replace('\n', '')
                    loc_names.append(line)

        while '' in loc_names:
            loc_names.remove('')

        sL = options['initState']
        startLoc = loc_names.index(sL) + 1
        res += "options.startLoc = " + str(startLoc) + ";\n"
        res += "options.finalLoc = 0;\n"
        res += "options.tStart = 0;\n"
        tFinal = options['tFinal']
        res += "options.tFinal = " + tFinal + ";\n"

        stepSize = options['stepSize']
        for i in range(1, len(loc_names) + 1):
            res += "options.timeStepLoc{" + str(i) + "}= " + stepSize + ";\n"

        res += "%define flows--------------------------------------------------------------\n\n"
        res += self.__constructFromConstraints(flows, options['stateVars'], 'flow', loc_names)
        res += "%define invariants---------------------------------------------------------\n\n"
        res += self.__constructFromConstraints(invariants, options['stateVars'], 'invariant', loc_names)

        inv_names = []
        return res, loc_names, inv_names

    def __constructFromConstraints(self, constraints, vars, name, loc_names):
        res = ""
        counter = 1
        for con in constraints:
            lhs = []
            operators = []
            rhs = []
            for c in con:
                c = c.replace('- ', '+ -')
                if '<' in c:
                    operators.append('<=')
                    c_array = c.split('<=')
                elif '>' in c:
                    operators.append('>=')
                    c_array = c.split('>=')
                elif '=' in c:
                    operators.append('=')
                    c_array = c.split('=')
                else:
                    sys.exit("Wrong operator in: " + c)

                lhs.append(c_array[0])
                rhs.append(c_array[1])

            if name == 'flow':
                normalized = self.__normalize(rhs, vars)
                A, b, intervals = self.__constraintsToMatrix(normalized, vars)
                A_matlab = self.__printMatrixToCORA(A)
                b_matlab = self.__printMatrixToCORA(b)
                res += 'A' + str(counter) + " = " + A_matlab + ';\n'
                res += 'B' + str(counter) + ' = zeros(' + str(len(vars)) + ", 1);\n"
                res += 'c' + str(counter) + " = " + b_matlab + ';\n'

                # process intervals
                inter_number = len([x for x in intervals if x != [0,0]])
                if inter_number > 0:
                    u_trans = []
                    delta = []
                    for i in intervals:
                        left = i[0]
                        right = i[1]
                        center = str((left + right) / 2)
                        d = str(abs(left - right))
                        u_trans.append(center)
                        delta.append(d)

                    res += "options.uTrans = " + self.__printVectorToCORA(u_trans) + ';\n'
                    res += "options.U = 0.5 * zonotope([zeros(" + str(len(lhs)) + ", 1), diag(" + self.__printVectorToCORA(delta) + ")]);\n"

                res += "flow" + str(counter) + " = linearSys('linearSys" + str(counter) + "', A" + str(
                    counter) + ", B" + str(counter) + ", c" + str(counter) + ");\n\n"
                counter += 1
            elif name == 'invariant':
                normalized = self.__normalize(lhs, vars)
                A, _, _ = self.__constraintsToMatrix(normalized, vars)
                res += 'inv_' + loc_names[counter] + " =  " + self.__writeMptPolytope(A, rhs, operators, vars)

            elif name == 'guard':
                normalized = self.__normalize(lhs, vars)
                A, _, _ = self.__constraintsToMatrix(normalized, vars)
                res += 'guard_' + loc_names[counter] + " =  " + self.__writeMptPolytope(A, rhs, operators, vars)

        return res

    def __normalize(self, expressions, vars):
        """
        This method normalizes the flow equations. For every missing component it adds a dummy.
        :rtype: object
        """
        normalized = []
        for expr in expressions:
            expr = expr.split(' + ')
            vars_contained = []
            constants = 0
            # Check which variables are contained in the expression
            for term in expr:
                term = term.replace(' ','')
                # Check if a constant is contained in the expression
                if self.__isNumber(term):
                    constants += 1
                else:
                    for v in vars:
                        if self.__isVariableInTerm(v,term):
                            vars_contained.append(v)

            #  Normalize
            if len(vars_contained) != len(vars):
                # Not all variables are contained in the expression
                for v in vars:
                    if v not in vars_contained:
                        # Add dummy
                        newTerm = 'non ' + v
                        expr.append(newTerm)
            if constants == 0:
                # There is no constant
                expr.append('0')
            if constants > 1:
                err_mes = "The expression for flow is not minimal! - Add the constants"
                sys.exit(err_mes)
            normalized.append(expr)
        return normalized

    def __constraintsToMatrix(self, constraints, vars):
        """
        This method converts a linear constraints into a Matlab format of matrix
        :param constraints: string containing the constraints
        :param vars: state variables
        :return: string containing the constraints in matrix form
        """
        matrix = []
        b = []
        intervals = []
        for f in constraints:
            entry = []
            coefficients = []
            # Get the coefficients for variables
            for v in vars:
                temp = [x for x in f if self.__isVariableInTerm(v,x)]
                if len(temp) > 1:
                    err_mes = 'The expression for constraints is not minimal! - Add the coefficients for variable ' + v
                    sys.exit(err_mes)
                if len(temp) == 0:
                    err_mes = "A coefficient in constraints is wrongly written!"
                    sys.exit(err_mes)
                c = temp[0].replace(' ', '')
                if ("-" + v) in c:
                    coefficients.append('-1')
                elif 'non' in c:
                    coefficients.append('0')
                elif '*' in c:
                    c_array = c.split('*')
                    coefficients.append(c_array[0])
                else:
                    coefficients.append('1')
            matrix.append(coefficients)

            # Get the constants
            temp = [x for x in f if self.__isNumber(x)]
            if len(temp) == 0:
                err_mes = "In one of the constraints expressions is no constant! - Add a 0"
                sys.exit(err_mes)
            b.append(temp)

            # Get intervals
            temp = [term for term in f if '[' in term]
            if len(temp) > 1:
                err_mes = "One constraints expression is wrong! - Added two intervals"
                sys.exit(err_mes)
            if len(temp) > 0:
                interval = temp[0]
                interval = interval.replace('[','')
                interval = interval.replace(']','')
                inter_array = interval.split(',')
                interval = [float(num_str) for num_str in inter_array]
                intervals.append(interval)
            elif len(temp) == 0:
                intervals.append([0,0])

        return matrix, b, intervals

    def __isVariableInTerm(self, var,term):
        term = term.replace(' ','')
        if '-' in term:
            term = term.replace('-','')
            term = term.replace(' ', '')
        if var == term:
            return True
        elif '*' in term:
            term_array = term.split('*')
            v = term_array[1]
            if var == v:
                return True
        elif 'non' in term:
            v = term.replace('non','')
            if var == v:
                return True
        else:
            return False


    def __writeMptPolytope(self, A, b, operators, vars):
        add_constr = []
        add_b = []

        print("A: ", A, " b: ", b)

        for i in range(len(operators)):
            if operators[i] == '>=':
                if b[i] != '0':
                    b[i] = '-' + b[i]
                for j in range(len(A[i])):
                    if A[i][j] != '0':
                        A[i][j] = '-' + A[i][j]
            elif operators[i] == '=':
                if b[i] != '0':
                    new_b = '-' + b[i]
                else:
                    new_b = '0'
                new_A_line = []
                for j in range(len(A[i])):
                    if A[i][j] != '0':
                        new_A_line.append('-' + A[i][j])
                    else:
                        new_A_line.append('0')
                add_constr.append(new_A_line)
                add_b.append(new_b)

        for item in add_constr:
            A.append(item)
        for item in add_b:
            b.append(item)

        A_matlab = self.__printMatrixToCORA(A)
        b_matlab = self.__printVectorToCORA(b)

        res = "mptPolytope(struct('A', "+ A_matlab + ", 'b', " + b_matlab + "));\n"

        return res

    def __printMatrixToCORA(self, matrix):
        res = "["
        cntr = 0
        print("matrix: ",matrix)
        print("rows: ", len(matrix))
        print("cols: ", len(matrix[0]))
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                res += str(matrix[i][j]) + " "
            cntr += 1
            if cntr < len(matrix):
                if len(matrix) > 5:
                    res += "; ...\n"
                else:
                    res += "; "
        res += "]"
        return res

    def __printVectorToCORA(self, vector):
        res = "["
        for v in range(len(vector)):
            res += vector[v]
            if v < len(vector) - 1:
                res += "; "
        res += "]"
        return res

    def __isNumber(self, expr):
        try:
            float(expr)
            return True
        except ValueError:
            return False

    def __getJumps(self, start, infile, vars, locations):
        res = "\n%define transitions--------------------------------------------------------\n"
        loc_dict = {}
        for l in locations:
            loc_dict[l] = 0

        open_braces = 1

        with open(infile, 'r') as file:
            for i in range(start + 1):
                file.readline()
            counter = 0
            l1 = ''
            l2 = ''
            for line in file:
                if open_braces == 0:
                    break
                if '{' in line:
                    open_braces += 1
                if '}' in line:
                    open_braces -= 1
                if '->' in line:
                    l = line.split('->')
                    l1 = l[0].replace(' ', '')
                    l1 = l1.replace('\n','')
                    l2 = l[1].replace(' ', '')
                    l2 = l2.replace('\n', '')
                    loc_dict[l1] += 1
                    res += '\n%' + line
                    counter += 1
                elif 'guard' in line:
                    line = line.replace('\n','')
                    line = line.replace('}','')
                    line = line.replace('{', '')
                    line = line.replace('guard', '')
                    guard_number = line.count('=')
                    if guard_number > 1:
                        line_array = line.split('   ')
                        while '' in line_array:
                            line_array.remove('')
                    else:
                        line_array = [line.replace(' ', '')]
                    self.__constructFromConstraints([line_array], vars, 'guard', locations)
                elif 'reset' in line:
                    if '{ }' in line:
                        res += "reset.A = eye(" + str(len(vars)) + ");\n"
                        res += "reset.b = zeros(" + str(len(vars)) + ", 1);\n"
                    else:
                        res += self.__extractConstraints(line, 'reset', vars)

                    counter += 1
                elif 'parallelotope' in line:
                    # TODO What is this?
                    pass
                else:
                    res += "trans_" + l1 + "{" + str(loc_dict[l1]) + "} = transition(guard" + ", reset, " +str(locations.index(l2) + 1) + ", '" + l1 + "', '" + l2 + "');\n"
        return res

    def __extractConstraints(self, line, name, vars, c = 0):
        res = ""
        line = line.replace('{','')
        line = line.replace('}', '')
        line = line.replace(name, '')
        line = line.replace('\n', '')
        line = line.replace('\t','')
        line = line.replace("'",'')
        line_array = line.split(' ')
        res += "%" + line +'\n'
        line_array = [x for x in line_array if len(x) > 0]

        counter = 0
        number_guards = 0
        lhs = []
        rhs = []
        operator = []
        for term in line_array:
            if counter == 0:
                if term in vars:
                    lhs.append(term)
                    number_guards += 1
                    counter += 1
                else:
                    sys.exit(" Undefined variable: " + term)
            elif counter == 1:
                if '=' in term or '>' in term or '<' in term:
                    operator.append(term)
                    counter +=1
                elif 'in' in term:
                    lhs.append(lhs[-1])
                    number_guards += 1
                    counter += 1
                    operator.append('>=')
                    operator.append('<=')
                else:
                    sys.exit("Constraint has to have the form [var] [operator] [constant] (wrong position of the operator): " + term)
            elif counter == 2:
                if self.__isNumber(term) or '*' in term or term in vars:
                    rhs.append(term)
                    counter = 0
                elif '[' in term:
                    term = term.replace('[','')
                    term = term.replace(']','')
                    interval = term.split(',')
                    rhs.append(interval[0])
                    rhs.append(interval[1])
                else:
                    sys.exit("Constraint has to have the form [var] [operator] [constant]: " +term)
            else:
                pass


        if name == 'guard' and number_guards == 0:
            #TODO is it possible in CORA?
            return res
        elif name == 'guard':
            res += "guard = " +  self.__writeMptPolytope(lhs, operator, rhs, vars)
        elif name == 'reset':
            res += self.__writeReset(lhs, operator, rhs,vars, c=counter)
        else:
            sys.exit("Wrong name")

        return res

    def __writeReset(self, lhs, operator, rhs, vars, c =0):
        res = ""
        A = []
        b = []
        var_indeces =[]
        for i in range(len(lhs)):
            if '>' in operator[i] or '<' in operator[i]:
                sys.exit("Reset have to be an :=")
            else:
                var_index = vars.index(lhs[i])
                var_indeces.append(var_index)
                for k in range(len(vars)):
                    entry = []
                    for v in range(len(vars)):
                        if v == var_index and k == var_index:
                            entry.append('0')
                        elif v == k :
                            entry.append('1')
                        else:
                            entry.append('0')
                    A.append(entry)

        counter = 0
        for i in range(len(vars)):
            if i in var_indeces:
                if self.__isNumber(rhs[counter]):
                    b.append(rhs[counter])
                elif '*' in rhs[counter]:
                    l = rhs[counter].split('*')
                    con = l[0].replace(' ', '')
                    A[i][i] = con
                    b.append('0')
                else:
                    A[i][i] = 1
                    b.append('0')
                counter += 1
            else:
                b.append('0')
        res += "reset.A = " + self.__printMatrixToCORA(A) + ";\n"
        res += "reset.b = " + self.__printVectorToCORA(b) + ";\n"
        return  res

    def __defineLocations(self, loc_names, inv_names):
        res = "\n%define locations----------------------------------------------------------\n\n"
        counter = 1
        for loc in loc_names:
            res += "options.uLoc{" + str(counter) + "} = 0;\n"
            res += "options.uLocTrans{" + str(counter) + "} = 0;\n"
            res += "options.Uloc{" + str(counter) + "} = zonotope(0);\n\n"
            counter += 1

        res += '\n'
        counter = 1

        for loc in loc_names:
            res += "loc{" + str(counter) + "} = location('" + loc_names[counter - 1] + "', " + str (counter) + ", " + "inv_" + loc + ", trans_"+ loc_names[counter - 1] + ", flow" + str(counter) + ");\n"
            counter += 1

        return res

    def __defineHybridAutomaton(self,reach=True, simulation=False):
        res = "\n%define hybrid automaton---------------------------------------------------\n\n"
        res += "HA = hybridAutomaton(loc);\n"
        res += "[HA] = reach(HA, options);\n\n"
        return res

    def __drawReachableSet(self):
        res = "figure\n"
        res += "hold on\n"
        res += "options.projectedDimensions = [1 2];\n"
        res += "options.plotType = 'b';\n"
        res += "plot(HA,'reachableSet',options); %plot reachable set\n"
        res += "plotFilled(options.R0,options.projectedDimensions,'w','EdgeColor','k'); %plot initial set\n\n"
        res += "res = 1;\n"
        return res


    def __writeOptions(self, options):
        res = ""
        res += "options.taylorTerms = " + str(options['taylor']) + ";\n"
        res += "options.zonotopeOrder = " + str(options['zonotope']) + ";\n"
        res += "options.polytopeOrder = " + str(options['polytope']) + ";\n"
        res += "options.reductionTechnique = '" + str(options['reduction']) + "';\n"
        res += "options.isHyperplaneMap = " + str(options['hyperplane']) + ";\n"
        res += "options.guardIntersect = '" + str(options['guard']) + "';\n"
        res += "options.enclosureEnables = " + str(options['enclosure']) + ";\n"
        res += "options.originContained = " + str(options['origin']) + ";\n"

        return res

    def convert(self, infile, outfile, options):
        # TODO check if which system we have
        path = infile.split('/')
        name = str(path[len(path) - 1])
        name = name.replace('.model', '')
        res = self.__readFile(infile, options, name)
        print(res)




class NonLinHybridFlowToCORA:
    pass


class LinContFlowToCORA:
    pass


class NonLinContFlowToCORA:
    pass
