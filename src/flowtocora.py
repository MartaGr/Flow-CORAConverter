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
        res += self.__getModes(modes_line, infile, options)
        res += self.__getJumps(jumps_line, infile)
        res += self.__getUnsafeSet(unsafe_line, infile)

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
                        line = line.split(' ')
                        it = line[len(line) - 1].replace('\n', '')
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
        res += self.__constructInvariants(loc_names, invariants, options['stateVars'])
        return res

    def __constructFlows(self,flows, vars):
        res = ""
        counter = 0

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
            interval = temp[0]

        m = self.__printMatrixToCORA(matrix)
        constants = self.__printMatrixToCORA(b)
        return m, constants, interval

    def __constructInvariants(self, loc_names, invariants, vars):
        #TODO Redo the invariants
        print("Invariants: ", invariants)


        return " "
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

    def __isNumber(self, expr):
        try:
            float(expr)
            return True
        except ValueError:
            return False

    def __getJumps(self, line, infile):
        res = ""
        return res

    def __getUnsafeSet(self, line, infile):
        res = ""
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
