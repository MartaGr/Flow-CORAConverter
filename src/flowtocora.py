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
                        line = line.split('=')
                        f = line[len(line) - 1].replace('\n', '')
                        fl.append(f)
                        line = file.readline()
                    flows.append(fl)
                if '{' in line and open_braces == 3 and inv_area:
                    # Save invariant
                    inv = []
                    line = file.readline()
                    while '}' not in line:
                        it = line.replace('\n', '')
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
        startLoc = loc_names.index(sL)
        res += "options.startLoc = " + str(startLoc) + ";\n"
        res += "options.finalLoc = 0;\n"
        res += "options.tStart = 0;\n"
        tFinal = options['tFinal']
        res += "options.tFinal = " + tFinal + ";\n"

        stepSize = options['stepSize']
        for i in range(1, len(loc_names)):
            res += "options.timeStepLoc{" + str(i) + "}= " + stepSize + ";\n\n"

        res += "%define flows--------------------------------------------------------------\n\n"
        res += self.__constructFlows(flows, options['stateVars'])
        res += "%define invariants---------------------------------------------------------\n\n"
        temp, inv_names= self.__constructInvariants(loc_names, invariants, options['stateVars'])
        res += temp
        return res, loc_names, inv_names

    def __constructFlows(self,flows, vars):
        res = ""
        counter = 1

        for flow in flows:
            single_flow = []
            for direc in flow:
                direc = direc.replace('- ', '+ -')
                single_flow.append(direc.split(' + '))
            normalized_flow = self.__normalize(single_flow, vars)
            A, b, inter = self.__flowToMatrix(normalized_flow, vars)
            res += 'A' + str(counter) + " = " + A + ';\n'
            res += 'B' + str(counter) + ' = zeros(' + str(len(vars)) + ", 1);\n"
            res += 'c' + str(counter) + " = " + b + ';\n'
            # TODO What to do with the intervals???
            res += "flow" + str(counter) + " = linearSys('linearSys" + str(counter) + "', A" + str(
                counter) + ", B" + str(counter) + ", c" + str(counter) + ");\n\n"
            counter += 1

        return res

    def __normalize(self, flow, vars):
        """
        This method normalizes the flow equations. For every missing component it adds a dummy.
        :rtype: object
        """
        for f in flow:
            vars_contained = []
            constants = 0
            # Check which variables are contained in the expression
            for term in f:
                for v in vars:
                    if v in term:
                        vars_contained.append(v)
                # Check if a constant is contained in the expression
                if self.__isNumber(term):
                    constants += 1

            #  Normalize
            if len(vars_contained) != len(vars):
                # Not all variables are contained in the expression
                for v in vars:
                    if v not in vars_contained:
                        # Add dummy
                        newTerm = '0 * ' + v
                        f.append(newTerm)
            if constants == 0:
                # There is no constant
                f.append('0')
            if constants > 1:
                errMes = "The expression for flow is not minimal! - Add the constants"
                sys.exit(errMes)
        return flow

    def __flowToMatrix(self, flow, vars):
        """
        This method converts a linear flow into a Matlab format of matrix
        :param flow: string containing the flow
        :param vars: state variables
        :return: string containing the flow in matrix form
        """
        matrix = []
        b = []
        interval = ''
        for f in flow:
            entry = []
            coefficients = []
            # Get the coefficients for variables
            for v in vars:
                temp = [x for x in f if (v in x)]
                if len(temp) > 1:
                    err_mes = 'The expression for flow is not minimal! - Add the coefficients for variable ' + v
                    sys.exit(err_mes)
                if len(temp) == 0:
                    err_mes = "A coefficient in flow is wrongly written!"
                    sys.exit(err_mes)
                c = temp[0].replace(' ', '')
                if ("-" + v) in c:
                    coefficients.append('-1')
                elif '0' in c:
                    coefficients.append('0')
                elif '*' in c:
                    c_array = c.split('*')
                    coefficients.append(c_array[0])
                else:
                    coefficients.append('1')
            # entry.append(coefficients)
            matrix.append(coefficients)

            # Get the constants
            temp = [x for x in f if self.__isNumber(x)]
            if len(temp) == 0:
                err_mes = "In one of the flow expressions is no constant! - Add a 0"
                sys.exit(err_mes)
            b.append(temp)

            # Get intervals
            temp = [term for term in f if '[' in term]
            if len(b) == 0:
                err_mes = "One flow expression is wrong! - Added two intervals"
                sys.exit(err_mes)
            if len(temp) > 0:
                interval = temp[0]

        m = self.__printMatrixToCORA(matrix)
        constants = self.__printMatrixToCORA(b)
        return m, constants, interval

    def __constructInvariants(self, loc_names, invariants, vars):
        res = ""
        counter = 0
        inv_names = []

        for inv in invariants:
            if inv != []:
                name = 'inv' + str(counter + 1)
                inv_names.append(name)
            else:
                inv_names.append('non')
            counter += 1
            lhs = []
            operator = []
            rhs = []
            for expr in inv:
                res += '%' + expr + '\n'
                i = expr.split(' ')
                if len(i) < 3:
                    err_mes = "An invariant has wrong form!"
                    sys.exit(err_mes)
                left = i[len(i) - 3]
                left = left.replace(' ', '')
                lhs.append(left)
                op = i[len(i) - 2]
                op = op.replace(' ', '')
                operator.append(op)
                right = i[len(i) - 1]
                right = right.replace(' ', '')
                rhs.append(right)
            print("invariant")
            res += 'inv' + str(counter) + " = " + self.__writeMptPolytope(lhs,operator,rhs, vars) + '\n'


        return res, inv_names

    def __writeMptPolytope(self, lhs, operator, rhs, vars):
        res = "mptPolytope(struct('A', "
        A = []
        b = []

        for i in range(len(lhs)):
            left = lhs[i].replace('\t','')
            left = left.replace('\n','')
            op = operator[i].replace('\t','')
            op = op.replace('\n', '')
            right = rhs[i].replace('\t','')
            right = right.replace('\n', '')
            print("left: ", left, " op:", op, " right: ", right)

            try:
                var_index = vars.index(left)
                constant = float(right)
            except ValueError:
                expr = left + op + right
                sys.exit("Wrong format: " + expr)
            if op == '=':
                b.append(str(-constant))
                b.append(str(constant))
                for k in range(2):
                    entry = []
                    for v in range(len(vars)):
                        if v == var_index and k == 0:
                            entry.append('-1')
                        elif v == var_index and k == 1:
                            entry.append('1')
                        else:
                            entry.append('0')
                    A.append(entry)
            elif op == "<=":
                b.append(str(constant))
                entry = []
                for v in range(len(vars)):
                    if v == var_index:
                        entry.append('1')
                    else:
                        entry.append('0')
                A.append(entry)
            elif op == ">=":
                b.append(str(-constant))
                entry = []
                for v in range(len(vars)):
                    if v == var_index:
                        entry.append('-1')
                    else:
                        entry.append('0')
                A.append(entry)
            else:
                err_mes = "One operand for invarinant is written wrong!"
                sys.exit(err_mes)

        if len(A) > 0:
            A_matlab = self.__printMatrixToCORA(A)
            b_matlab = self.__printVectorToCORA(b)
            print(A_matlab + ", 'b', " + b_matlab)
            res += A_matlab + ", 'b', " + b_matlab + "));\n"
        else:
            # There is no guard or invariant
            # TODO what does CORA do with this???
            res += ""
        return res

    def __printMatrixToCORA(self, matrix):
        res = "["
        cntr = 0
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
                    res += self.__extractGuards(line, vars)
                elif 'reset' in line:
                    if '{}' in line:
                        res += "reset" + str(counter) + ".A = eye(" + str(len(vars)) + ");\n"
                        res += "reset" + str(counter) + ".b = zeros(" + str(len(vars)) + ", 1);\n"
                elif 'parallelotope' in line:
                    # TODO What is this?
                    pass
                else:
                    res += "trans_" + l1 + "{" + str(loc_dict[l1]) + "} = transition(guard" + str(counter) + ", reset" + str(counter) + ", " +str(locations.index(l2) + 1) + ", '" + l1 + "', '" + l2 + "');\n"
        return res

    def __extractGuards(self, line, vars):
        res = ""
        line = line.replace('{','')
        line = line.replace('}', '')
        line = line.replace('guard', '')
        line = line.replace('\n', '')
        line = line.replace('\t','')
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
                    sys.exit("Guard error: Undefined variable: " + term)
            elif counter == 1:
                if '=' in term or '>' in term or '<' in term:
                    operator.append(term)
                    counter +=1
                else:
                    sys.exit("Guard error: the guard has to have the form [var] [operator] [constant]")
            elif counter == 2:
                if self.__isNumber(term):
                    rhs.append(term)
                    counter = 0
                else:
                    sys.exit("Guard error: the guard has to have the form [var] [operator] [constant]")
            else:
                pass


        if number_guards == 0:
            #TODO is it possible in CORA?
            return res
        else:
            res += "guard = " +  self.__writeMptPolytope(lhs, operator, rhs, vars)

        return res

    def __defineLocations(self, loc_names, inv_names):
        res = "\n%define locations----------------------------------------------------------\n\n"
        counter = 1
        for loc in loc_names:
            res += "options.uLoc{" + str(counter) + "} = 0;\n"
            res += "options.uLocTrans{" + str(counter) + "} = 0;\n"
            res += "options.Uloc{" + str(counter) + "} = 0;\n\n"
            counter += 1

        res += '\n'
        counter = 1

        for loc in loc_names:
            if inv_names[counter - 1] != "non":
                inv = "inv" + str(counter)
            else:
                inv = '0'
            res += "loc{" + str(counter) + "} = location('" + loc_names[counter - 1] + "', " + str (counter) + ", " + inv + ", trans_"+ loc_names[counter - 1] + ", flow" + str(counter) + ");\n"
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
        res += "options.reductionTechnique = " + str(options['reduction']) + ";\n"
        res += "options.isHyperplaneMap = " + str(options['hyperplane']) + ";\n"
        res += "options.guardIntersect = " + str(options['guard']) + ";\n"
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
