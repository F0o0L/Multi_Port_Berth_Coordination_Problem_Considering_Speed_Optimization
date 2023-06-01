import cplex
import numpy as np
import time
from Initialization import Initialization


def CplexBAP(ini: Initialization):
    start = time.time()
    M = 10000

    T = []
    a = []
    deltaEFT = []
    t = []
    q = []
    u = []
    Fc = []
    tlb = []
    tub = []

    IcNum1 = 0
    IcNum2 = 0
    for ves, vessel in enumerate(ini.portOrder):
        for p, port in enumerate(vessel):
            IcNum2 += 1
            a.append('a:i{}_p{}_port{}'.format(ves, p, port))
            deltaEFT.append('deltaEFT:i{}_p{}_port{}'.format(ves, p, port))
            t.append('t:i{}_p{}_port{}'.format(ves, p, port))
            q.append('q:i{}_p{}_port{}'.format(ves, p, port))
            u.append('u:i{}_p{}_port{}_num{}'.format(ves, p, port, 1))
            u.append('u:i{}_p{}_port{}_num{}'.format(ves, p, port, 2))
            if p == 0:
                Fc.append(ini.Fc * ini.fd * np.power(ini.distances[1][ves] / ini.speedLimit[0], 3))
            else:
                Fc.append(ini.Fc * ini.fd * np.power(ini.distances[0][vessel[p - 1], port] / ini.speedLimit[0], 3))
            tlb.append(ini.travelTimeLimit[ves][p, 0])
            tub.append(ini.travelTimeLimit[ves][p, 1])

            for berth in range(ini.berthNum[port]):
                IcNum1 += 1
                T.append('T:ves{}_p{}_port{}_berth{}'.format(ves, p, port, berth))

    Ic1 = [ini.Ic] * IcNum1
    Ic2 = [-ini.Ic] * IcNum2
    Dc = [ini.Dc] * IcNum2

    x = []
    Hc = []
    for p in range(ini.portNum):
        vessels = ini.vesselsInPort[p]
        for v1 in vessels:
            for v2 in vessels:
                if v1 != v2:
                    for berth in range(ini.berthNum[p]):
                        Hc.append(ini.handlingTimes[p][berth, v1] * ini.Hc)
                        x.append('x:i{}_j{}_port{}_berth{}'.format(v1, v2, p, berth))

    IcNum1 = 0
    for p, port in enumerate(ini.vesselsInPort):
        for berth in range(ini.berthNum[p]):
            for ves in port:
                x.append('x:i{}_j{}'.format('o(port' + str(p) + ',berth' + str(berth) + ')', ves))
                Hc.append(0)
                x.append('x:i{}_j{}'.format(ves, 'd(port' + str(p) + ',berth' + str(berth) + ')'))
                Hc.append(ini.handlingTimes[p][berth, ves] * ini.Hc)
            x.append('x:i{}_j{}'.format('o(port' + str(p) + ',berth' + str(berth) + ')',
                                        'd(port' + str(p) + ',berth' + str(berth) + ')'))
            Hc.append(0)
            T.append('T:o(port' + str(p) + ',berth' + str(berth) + ')')
            T.append('T:d(port' + str(p) + ',berth' + str(berth) + ')')
            IcNum1 += 1
    Ic1 += [0] * IcNum1 * 2

    obj = Ic1 + Ic2 + Hc + Dc + Fc + [0] * (len(u + t))
    lb = [0] * len(Ic1 + Ic2 + Hc + Dc + Fc + u) + tlb
    ub = [cplex.infinity] * (len(Ic1 + Ic2)) + [1] * len(Hc) + [cplex.infinity] * (len(Dc + Fc + u)) + tub
    cType = 'C' * (len(Ic1 + Ic2)) + 'I' * len(Hc) + 'C' * (len(Dc + Fc + u + t))
    colNames = T + a + x + deltaEFT + q + u + t

    rows = []
    rowNames = []
    rhs = []
    sense = ''
    socpRows = []
    socpRowNames = []
    socpRhs = []
    socpSense = ''

    for i in range(ini.vesselNum):
        for p, port in enumerate(ini.portOrder[i]):
            row = []
            for st in x:
                if 'x:i' + str(i) + '_' in st:
                    if '_port' + str(port) in st:
                        row.append(st)
                    if 'jd(port' + str(port) in st:
                        row.append(st)
            rows.append([row, [1] * len(row)])
            rowNames.append('1_i' + str(i) + '_p' + str(p) + '_port' + str(port))
            rhs.append(1)
            sense += 'E'

            for berth in range(ini.berthNum[port]):
                row1 = []
                row2 = []
                num1 = []
                num2 = []
                row2.append('T:ves{}_p{}_port{}_berth{}'.format(i, p, port, berth))
                num2.append(1)
                for st in x:
                    if 'x:i' + str(i) + '_' in st:
                        if '_port' + str(port) in st:
                            if '_berth' + str(berth) in st:
                                t1 = int(st[st.find('x:i') + 3:st.find('_j')])
                                t2 = int(st[st.find('_j') + 2:st.find('_port')])
                                if t1 != t2:
                                    row1.append(st)
                                    num1.append(1)
                                row2.append(st)
                                num2.append(-M)
                        if 'jd(port' + str(port) + ',berth' + str(berth) + ')' in st:
                            row1.append(st)
                            num1.append(1)
                            row2.append(st)
                            num2.append(-M)
                    tempIdx = st.find('_j' + str(i))
                    if '_j' + str(i) + '_' in st or (tempIdx != -1 and tempIdx + len('_j' + str(i)) == len(st)):
                        if '_port' + str(port) in st:
                            if '_berth' + str(berth) in st:
                                t1 = int(st[st.find('x:i') + 3:st.find('_j')])
                                t2 = int(st[st.find('_j') + 2:st.find('_port')])
                                if t1 != t2:
                                    row1.append(st)
                                    num1.append(-1)
                                    row2.append(st)
                                    num2.append(-M)
                        if 'io(port' + str(port) + ',berth' + str(berth) + ')' in st:
                            row1.append(st)
                            num1.append(-1)
                            row2.append(st)
                            num2.append(-M)
                rows.append([row1, num1])
                rowNames.append('4_i' + str(i) + '_p' + str(p) + '_port' + str(port) + '_berth' + str(berth))
                rhs.append(0)
                sense += 'E'
                rows.append([row2, num2])
                rowNames.append('11_i' + str(i) + '_p' + str(p) + '_port' + str(port) + '_berth' + str(berth))
                rhs.append(0)
                sense += 'L'

            if p == 0:
                row = []
                num = []
                row.append('t:i{}_p{}_port{}'.format(i, p, port))
                num.append(1)
                row.append('a:i{}_p{}_port{}'.format(i, p, port))
                num.append(-1)
                rows.append([row, num])
                rowNames.append('6i' + str(i) + '_p' + str(p) + '_port' + str(port))
                rhs.append(0)
                sense += 'E'
            if len(ini.portOrder[i]) - 1 > p:
                row = []
                num = []
                for st in T:
                    if 'T:ves{}_p{}_port{}'.format(i, p, port) in st:
                        row.append(st)
                        num.append(1)
                for st in x:
                    if 'x:i{}'.format(i) + '_' in st:
                        if '_port' + str(port) in st or '_j{}'.format('d(port' + str(p)) in st:
                            t1 = st.find('berth') + 5
                            t2 = st.find(')')
                            if t2 != -1:
                                berth = int(st[t1:t2])
                            else:
                                berth = int(st[t1:len(st)])
                            row.append(st)
                            num.append(ini.handlingTimes[port][berth, i])
                row.append('t:i{}_p{}_port{}'.format(i, p + 1, ini.portOrder[i][p + 1]))
                num.append(1)
                row.append('a:i{}_p{}_port{}'.format(i, p + 1, ini.portOrder[i][p + 1]))
                num.append(-1)
                rows.append([row, num])
                rowNames.append(
                    '7i' + str(i) + '_p' + str(p) + '_port' + str(port) + '_pNext' + str(p + 1) + '_portNext' + str(
                        ini.portOrder[i][p + 1]))
                rhs.append(0)
                sense += 'E'

            row = []
            num = []
            for st in T:
                if 'T:ves{}_p{}_port{}'.format(i, p, port) in st:
                    row.append(st)
                    num.append(1)
            rows.append([row, num])
            rowNames.append('8i' + str(i) + '_p' + str(p) + '_port' + str(port))
            rhs.append(ini.timeWindows[i][p, 0])
            sense += 'G'

            row1 = []
            num1 = []
            row1.append('deltaEFT:i{}_p{}_port{}'.format(i, p, port))
            num1.append(1)
            row1.extend(row)
            num1.extend([-1] * len(row))

            row2 = []
            num2 = []
            row2.append('a:i{}_p{}_port{}'.format(i, p, port))
            num2.append(1)
            row2.extend(row)
            num2.extend([-1] * len(row))
            for st in x:
                if 'x:i{}'.format(i) + '_' in st:
                    if '_port' + str(port) in st or '_j{}'.format('d(port' + str(p)) in st:
                        t1 = st.find('berth') + 5
                        t2 = st.find(')')
                        if t2 != -1:
                            berth = int(st[t1:t2])
                        else:
                            berth = int(st[t1:len(st)])
                        row1.append(st)
                        num1.append(-ini.handlingTimes[port][berth, i])
            rows.append([row1, num1])
            rowNames.append('9i' + str(i) + '_p' + str(p) + '_port' + str(port))
            rhs.append(-ini.timeWindows[i][p, 1])
            sense += 'G'
            rows.append([row2, num2])
            rowNames.append('10i' + str(i) + '_p' + str(p) + '_port' + str(port))
            rhs.append(0)
            sense += 'L'

            socpRows.append([[['u:i{}_p{}_port{}_num{}'.format(i, p, port, 1)],
                              ['u:i{}_p{}_port{}_num{}'.format(i, p, port, 1)], [1]],
                             [['q:i{}_p{}_port{}'.format(i, p, port)], [-1]]])
            socpRowNames.append('socp_i{}_p{}_port{}_num{}'.format(i, p, port, 1))
            socpRhs.append(0)
            socpSense += 'L'

            socpRows.append([[['u:i{}_p{}_port{}_num{}'.format(i, p, port, 1)],
                              ['t:i{}_p{}_port{}'.format(i, p, port)], [1]], []])
            socpRowNames.append('socp_i{}_p{}_port{}_num{}'.format(i, p, port, 2))
            socpRhs.append(1)
            socpSense += 'G'

    for port in range(ini.portNum):
        for k in range(ini.berthNum[port]):
            row1 = []
            row2 = []
            for st in x:
                if 'io(port' + str(port) + ',berth' + str(k) + ')' in st:
                    row1.append(st)
                if 'jd(port' + str(port) + ',berth' + str(k) + ')' in st:
                    row2.append(st)
            rows.append([row1, [1] * len(row1)])
            rowNames.append('2_port' + str(port) + '_berth' + str(k))
            rhs.append(1)
            sense += 'E'
            rows.append([row2, [1] * len(row2)])
            rowNames.append('3_port' + str(port) + '_berth' + str(k))
            rhs.append(1)
            sense += 'E'

    for p, vessels in enumerate(ini.vesselsInPort):
        for k in range(ini.berthNum[p]):
            row = []
            num = []
            t1 = 'T:o(port' + str(p) + ',berth' + str(k) + ')'
            t2 = 'T:d(port' + str(p) + ',berth' + str(k) + ')'
            x1 = 'x:i{}_j{}'.format('o(port' + str(p) + ',berth' + str(k) + ')',
                                    'd(port' + str(p) + ',berth' + str(k) + ')')
            row.append(t1)
            row.append(t2)
            row.append(x1)
            num.append(1)
            num.append(-1)
            num.append(M)
            rows.append([row, num])
            rowNames.append(
                '5i' + 'o(port' + str(p) + ',berth' + str(k) + ')' + '_j' + 'd(port' + str(p) + ',berth' + str(k) + ')')
            rhs.append(M)
            sense += 'L'

            for i in vessels:
                row1 = []
                num1 = []
                t1 = 'T:o(port' + str(p) + ',berth' + str(k) + ')'
                t2 = 'T:ves{}_p{}_port{}_berth{}'.format(i, int(np.where(ini.portOrder[i] == p)[0][0]), p, k)
                x1 = 'x:i{}_j{}'.format('o(port' + str(p) + ',berth' + str(k) + ')', i)
                row1.append(t1)
                row1.append(t2)
                row1.append(x1)
                num1.append(1)
                num1.append(-1)
                num1.append(M)
                rows.append([row1, num1])
                rowNames.append('5i' + 'o(port' + str(p) + ',berth' + str(k) + ')' + '_j' + str(i))
                rhs.append(M)
                sense += 'L'

                row2 = []
                num2 = []
                t1 = 'T:ves{}_p{}_port{}_berth{}'.format(i, int(np.where(ini.portOrder[i] == p)[0][0]), p, k)
                t2 = 'T:d(port' + str(p) + ',berth' + str(k) + ')'
                x1 = 'x:i{}_j{}'.format(i, 'd(port' + str(p) + ',berth' + str(k) + ')')
                row2.append(t1)
                row2.append(t2)
                row2.append(x1)
                num2.append(1)
                num2.append(-1)
                num2.append(M)
                rows.append([row2, num2])
                rowNames.append('5i' + str(i) + '_j' + 'd(port' + str(p) + ',berth' + str(k) + ')')
                rhs.append(M - ini.handlingTimes[p][k, i])
                sense += 'L'

            for i in vessels:
                for j in vessels:
                    if i != j:
                        row = []
                        num = []
                        t1 = 'T:ves{}_p{}_port{}_berth{}'.format(i, int(np.where(ini.portOrder[i] == p)[0][0]), p, k)
                        t2 = 'T:ves{}_p{}_port{}_berth{}'.format(j, int(np.where(ini.portOrder[j] == p)[0][0]), p, k)
                        x1 = 'x:i{}_j{}_port{}_berth{}'.format(i, j, p, k)
                        row.append(t1)
                        row.append(t2)
                        row.append(x1)
                        num.append(1)
                        num.append(-1)
                        num.append(M)
                        rows.append([row, num])
                        rowNames.append('5i' + str(i) + '_j' + str(j) + '_port' + str(p) + '_berth' + str(k))
                        rhs.append(M - ini.handlingTimes[p][k, i])
                        sense += 'L'
    myProb = cplex.Cplex()
    myProb.parameters.timelimit.set(10800)
    myProb.objective.set_sense(myProb.objective.sense.minimize)
    myProb.variables.add(obj=obj, lb=lb, ub=ub, types=cType, names=colNames)
    myProb.linear_constraints.add(lin_expr=rows, senses=sense, rhs=rhs, names=rowNames)
    for so in range(len(socpRowNames)):
        if len(socpRows[so][1]) == 0:
            myProb.quadratic_constraints.add(quad_expr=socpRows[so][0], sense=socpSense[so], rhs=socpRhs[so],
                                             name=socpRowNames[so])
        else:
            myProb.quadratic_constraints.add(lin_expr=socpRows[so][1], quad_expr=socpRows[so][0],
                                             sense=socpSense[so],
                                             rhs=socpRhs[so],
                                             name=socpRowNames[so])
    myProb.write('MPBAP.lp')
    myProb.solve()
    end = time.time()
    runTime = end - start
    lb = myProb.solution.MIP.get_best_objective()
    gap = myProb.solution.MIP.get_mip_relative_gap()
    print(myProb.solution.get_objective_value())
    # values = myProb.solution.get_values()
    # for i in range(len(colNames)):
    #     print(colNames[i] + ' ' + cType[i] + ' ' + str(values[i]))

    berthAllocation = []
    for _ in range(ini.portNum):
        temp = []
        for __ in range(ini.berthNum[_]):
            temp.append([])
        berthAllocation.append(temp)

    xSolution = dict(zip(x, myProb.solution.get_values(x)))
    tSolution = dict(zip(t, myProb.solution.get_values(t)))
    TSolution = dict(zip(T, myProb.solution.get_values(T)))
    for p, port in enumerate(ini.vesselsInPort):
        for berth in range(ini.berthNum[p]):
            for ves in port:
                if int(np.round(
                        xSolution[
                            'x:i{}_j{}'.format('o(port' + str(p) + ',berth' + str(berth) + ')', ves)])) == 1:
                    berthAllocation[p][berth].append(ves)

    for p in range(ini.portNum):
        for b in range(ini.berthNum[p]):
            if len(berthAllocation[p][b]) != 0:
                temp = berthAllocation[p][b][len(berthAllocation[p][b]) - 1]
                while True:
                    for ves in range(ini.vesselNum):
                        if ves != temp:
                            if int(np.round(
                                    xSolution['x:i{}_j{}_port{}_berth{}'.format(temp, ves, p, b)])) == 1:
                                berthAllocation[p][b].append(ves)
                                temp = ves
                                break
                    else:
                        break
    travelTime = []
    for ves in range(ini.vesselNum):
        travelTime.append(np.zeros(len(ini.portOrder[ves])))

    for ves in range(ini.vesselNum):
        for p, port in enumerate(ini.portOrder[ves]):
            travelTime[ves][p] = tSolution['t:i{}_p{}_port{}'.format(ves, p, port)]

    handlingBeginTimeTable = []

    for ves in range(ini.vesselNum):
        handlingBeginTimeTable.append([])

    for ves in range(ini.vesselNum):
        for p, port in enumerate(ini.portOrder[ves]):
            for b in range(ini.berthNum[port]):
                temp = TSolution['T:ves{}_p{}_port{}_berth{}'.format(ves, p, port, b)]
                if temp != 0:
                    handlingBeginTimeTable[ves].append(temp)
                    break

    totalCost, travelTime, arriveTable, handlingBeginTimeTable, handlingEndTimeTable = ini.calculateObjectiveFunction(
        travelTime=travelTime,
        tryBerthAllocation=berthAllocation,
        # handlingBeginTimeTable=handlingBeginTimeTable,
        outputTime=True)
    print(totalCost)
    print(berthAllocation)
    print(travelTime)
    print(arriveTable)
    print(handlingBeginTimeTable)
    print(handlingEndTimeTable)
    # ini.paintSolution(travelTime, berthAllocation, handlingBeginTimeTable=handlingBeginTimeTable)

    with open('result.txt', 'a') as f:
        f.write('Cplex_portNum{}_berthNum{}_vesselNum{}\n'.format(ini.portNum, ini.berthNum, ini.vesselNum))
        f.write('totalCost:{}\n'.format(myProb.solution.get_objective_value()))
        f.write('lb:{}\n'.format(lb))
        f.write('gap:{}\n'.format(gap))
        f.write('berthAllocation:{}\n'.format(berthAllocation))
        travelTime = str(travelTime).replace('array', 'np.array')
        f.write('travelTime:{}\n'.format(travelTime))
        f.write('runTime:{}\n'.format(runTime))


if __name__ == '__main__':
    ini = Initialization(portNum=3, berthNum=3, vesselNum=5)
    CplexBAP(ini)
