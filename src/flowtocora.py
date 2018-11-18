import sys
import re

class LinHybridFlowStarToCORA:

    def __readFile(self, infile, options, name):
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
                    stepSize = stepSize.replace('\n','')
                    options['stepSize'] = stepSize
                if 'time' in line:
                    tFinal = line.split(' ')[3]
                    tFinal = tFinal.replace('\n','')
                    options['tFinal'] = tFinal
                if 'remainder estimation' in line:
                    reminder = line.split(' ')[4]
                    reminder = reminder.replace('\n','')
                    options['reminder'] = reminder
                if 'max jumps' in line:
                    jumps = line.split(' ')[4]
                    jumps = jumps.replace('\n','')
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
        res += self.__getInitialization(init_line,infile,name,options)
        res += self.__writeOptions(options)
        res += self.__getModes(modes_line,infile,options)
        res += self.__getJumps(jumps_line,infile)
        res += self.__getUnsafeSet(unsafe_line,infile)

        return res


    def __getInitialization(self, start, infile,name,options):
        res = "function res = " + name + "()\n\n"
        res += "%set options---------------------------------------------------------------\n"
        initial = []

        x0 = "options.x0 = ["
        r0 = "options.R0 = zonotope([options.x0, "

        with open(infile, 'r') as file:
            for i in range(start+2):
                file.readline()
                if i == start:
                    initState = file.readline()
                    initState = initState.replace('\n','')
                    initState = initState.replace(' ', '')
                    options['initState'] = initState

            #Now, extract the initalization information
            for line in file:
                if '}' in line:
                    break
                else:
                    l = line.split(' ')
                    initial.append(l[len(l) - 1])
        file.close()

        count = 0
        for item in initial:
            item = item.replace('[','')
            item = item.replace(']','')
            limits = item.split(',')
            lower = float(limits[0])
            upper = float(limits[1])
            middle = (upper - lower)/2
            x0 += str(lower + middle)
            if count < len(initial) - 1:
                x0 += "; "
            count += 1


        x0 += "];\n"
        r0 += options['stepSize'] + " * eye(" + str(len(initial)) + ")]);\n"

        res += x0
        res += r0
        return res

    def __getModes(self, start, infile,options):
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
                        f = line[len(line) - 1].replace('\n','')
                        fl.append(f)
                        line = file.readline()
                    flows.append(fl)
                if '{' in line and open_braces == 3 and inv_area:
                    # Save invariant
                    inv = []
                    line = file.readline()
                    while '}' not in line:
                        line = line.split(' ')
                        it = line[len(line) - 1].replace('\n','')
                        inv.append(it)
                        line = file.readline()
                    invariants.append(inv)
                    inv_area = False
                    open_braces -= 1
                if open_braces == 1 and '}' not in line:
                    line = line.replace(' ', '')
                    line = line.replace('\n','')
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
        for i in range(1,len(loc_names)):
            res += "options.timeStepLoc{" + str(i) + "}= " + stepSize + ";\n"

        res += "\n"
        res += self.__constructLocations(loc_names, flows, invariants, options)

        return res

    def __constructLocations(self, loc_names, flows, invariants, options):
        res = "%define flows--------------------------------------------------------------\n"
        vars = options['stateVars']
        intervals = []


        for flow in flows:
            print("Next flow")
            single_flows = []
            matrix = []
            b = []
            # for direc in flow:
            #     direc = direc.replace('- ', '+ -')
            #     single_flows.append(direc.split(' + '))
            # print(single_flows)
            for fl in single_flows:
                print("Fl: ", fl)
                for term in fl:
                    # Check if it is an interval
                    if '[' in term:
                        intervals.append(term)
                    else:
                        # Check if a variable is contained
                        contains_var = False
                        for v in vars:
                            if v in term:
                                # The term contains a variable
                                term = term.replace(v, '')
                                term = term.replace('*','')
                                term = term.replace(' ','')
                                if not any(char.isdigit() for char in term):
                                    term = '-1'
                                entry = []
                                var_index = vars.index(v)
                                for i in range(len(vars)):
                                    if i != var_index:
                                        entry.append(0)
                                    else:
                                        entry.append(term)
                                print("New entry: ",entry)
                                matrix.append(entry)
                                contains_var = True
                                break
                        if not contains_var:
                            # The term is just a constant
                            term = term.replace(' ', '')
                            entry = []
                            for i in range(len(vars)):
                                entry.append(0)
                            b.append(term)
                            matrix.append(entry)

            print('A: ', matrix)
            print('b: ', b)
        return res




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
        #TODO check if which system we have
        path = infile.split('/')
        name = str(path[len(path) - 1])
        name = name.replace('.model','')
        res = self.__readFile(infile, options, name)
        print(res)


class NonLinHybridFlowToCORA:
    pass

class LinContFlowToCORA:
    pass

class NonLinContFlowToCORA:
    pass

