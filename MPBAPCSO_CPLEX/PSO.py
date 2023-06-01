import numpy as np
import copy
import time
from Initialization import Initialization


class PSO:
    def __init__(self, ini: Initialization, swarmNum=200, iterationNum=2000):
        self.ini = ini

        self.swarmNum = swarmNum
        self.iterationNum = iterationNum
        self.c1 = 2
        self.c2 = 2
        self.sumPorts = 0
        for ports in self.ini.portOrder:
            self.sumPorts += len(ports)

        self.travelTimeVelocityLimit = self.getTravelTimeVelocityLimit()

        Pc = np.arange(1, self.swarmNum + 1)
        self.Pc = 0.05 + 0.45 * ((np.exp(10 * (Pc - 1) / (self.swarmNum - 1)) - 1) / (np.exp(10) - 1))
        self.rg = 7

        self.swarm = self.initSwarm()
        self.GBest = None
        self.PBest = None

        self.runtime = None

    def initSwarm(self):
        swarm = []
        for _ in range(self.swarmNum):
            solution = self.ini.generateFeasibleSolution()
            swarm.append((solution[0], solution[1], self.berthAllocationTransferToArray(solution[2])))
        return swarm

    def getGBest(self):
        idx = 0
        for i in range(self.swarmNum):
            if self.swarm[i][0] < self.swarm[idx][0]:
                idx = i
        if self.GBest is None:
            self.GBest = copy.deepcopy(self.swarm[idx])
        else:
            if self.swarm[idx][0] < self.GBest[0]:
                self.GBest = copy.deepcopy(self.swarm[idx])

    def getPBest(self):
        if self.PBest is None:
            self.PBest = copy.deepcopy(self.swarm)
        else:
            for i in range(self.swarmNum):
                if self.swarm[i][0] < self.PBest[i][0]:
                    self.PBest[i] = copy.deepcopy(self.swarm[i])

    def berthAllocationTransferToArray(self, berthAllocation):
        berthAllocationArray = []
        for bn, ba, v in zip(self.ini.berthNum, berthAllocation, self.ini.vesselsInPort):
            temp = np.zeros((2 * bn + len(v), 2 * bn + len(v)))
            for i, b in enumerate(ba):
                for j, ve in enumerate(b):
                    if j == 0:
                        temp[i * 2, 2 * bn + v.index(ve)] = 1
                    if j == len(b) - 1:
                        temp[2 * bn + v.index(ve), i * 2 + 1] = 1
                    elif j <= len(b) - 1:
                        temp[2 * bn + v.index(ve), 2 * bn + v.index(b[j + 1])] = 1
            berthAllocationArray.append(temp)
        return berthAllocationArray

    def berthsTransferToList(self, port, berthsArray):
        bn = self.ini.berthNum[port]
        v = self.ini.vesselsInPort[port]
        berths = []
        for i in range(bn):
            temp = []
            idx = np.where(berthsArray[2 * i, :] == 1)
            if len(idx[0]) > 0:
                idx = idx[0][0]
                while idx != 2 * i + 1:
                    temp.append(v[idx - 2 * bn])
                    idx = np.where(berthsArray[idx, :] == 1)[0][0]
            berths.append(temp)
        return berths

    def berthAllocationTransferToList(self, berthAllocationArray):
        berthAllocation = []
        for i in range(self.ini.portNum):
            berthAllocation.append(self.berthsTransferToList(i, berthAllocationArray[i]))
        return berthAllocation

    def getTravelTimeVelocityLimit(self):
        travelTimeVelocityLimit = []
        for ves in range(self.ini.vesselNum):
            travelTimeVelocityLimit.append([])

        for ves in range(self.ini.vesselNum):
            for p in range(len(self.ini.portOrder[ves])):
                limit = (self.ini.travelTimeLimit[ves][p, 1] - self.ini.travelTimeLimit[ves][p, 0]) * 0.15
                travelTimeVelocityLimit[ves].append(limit)

        return travelTimeVelocityLimit

    def initVelocity(self):
        travelTimeVelocity = []
        for ves in range(self.ini.vesselNum):
            travelTimeVelocity.append([])

        for ves in range(self.ini.vesselNum):
            for p in range(len(self.ini.portOrder[ves])):
                limit = self.travelTimeVelocityLimit[ves][p]
                travelTimeVelocity[ves].append(np.random.uniform(-limit, limit))

        berthAllocationVelocity = []
        for i, berthNum in enumerate(self.ini.berthNum):
            vesselNumInPort = len(self.ini.vesselsInPort[i])
            tableShape = 2 * berthNum + vesselNumInPort
            temp = np.zeros((tableShape, tableShape, 2))
            n = np.random.randint(vesselNumInPort + 1, vesselNumInPort * 2 - max(vesselNumInPort - berthNum, 0) + 1)
            for _ in range(n):
                r1List = []
                r2List = []
                for j in range(0, 2 * berthNum, 2):
                    r1List.append(j)
                    r2List.append(j + 1)
                for j in range(2 * berthNum, tableShape):
                    r1List.append(j)
                    r2List.append(j)
                r1 = np.random.choice(r1List)
                r2 = np.random.choice(r2List)
                while r1 < 2 * berthNum and r2 < 2 * berthNum or r1 == r2 or temp[r1, r2, 0] == 1:
                    r1 = np.random.choice(r1List)
                    r2 = np.random.choice(r2List)
                temp[r1, r2, 0] = 1
                temp[r1, r2, 1] = np.random.rand()
            berthAllocationVelocity.append(temp)
        return travelTimeVelocity, berthAllocationVelocity

    def CLPSOPositionSubtraction(self, c, idx):
        travelTimeVel = []
        for ves in range(self.ini.vesselNum):
            travelTimeVel.append([])

        for ves in range(self.ini.vesselNum):
            for p in range(len(self.ini.portOrder[ves])):
                Pc = np.random.rand()
                if Pc > self.Pc[idx]:
                    newVel = c * np.random.rand() * (self.PBest[idx][1][ves][p] - self.swarm[idx][1][ves][p])
                else:
                    idx1 = np.random.randint(0, self.swarmNum)
                    idx2 = np.random.randint(0, self.swarmNum)
                    if idx1 == idx or idx2 == idx or idx1 == idx2:
                        idx1 = np.random.randint(0, self.swarmNum)
                        idx2 = np.random.randint(0, self.swarmNum)
                    if self.swarm[idx1][0] < self.swarm[idx2][0]:
                        newVel = c * np.random.rand() * (self.PBest[idx1][1][ves][p] - self.swarm[idx][1][ves][p])
                    else:
                        newVel = c * np.random.rand() * (self.PBest[idx2][1][ves][p] - self.swarm[idx][1][ves][p])
                travelTimeVel[ves].append(newVel)

        berthAllocationVel = []
        for port in range(self.ini.portNum):
            num = self.ini.berthNum[port] * 2 + len(self.ini.vesselsInPort[port])
            berthAllocationVel.append(np.zeros((num, num, 2)))

        for port in range(self.ini.portNum):
            berthNum = self.ini.berthNum[port]
            for berth in range(self.ini.berthNum[port]):
                Pc = np.random.rand()
                if Pc > self.Pc[idx]:
                    newVel = self.PBest[idx][2][port][2 * berth, :] - self.swarm[idx][2][port][2 * berth, :]
                else:
                    idx1 = np.random.randint(0, self.swarmNum)
                    idx2 = np.random.randint(0, self.swarmNum)
                    if idx1 == idx or idx2 == idx or idx1 == idx2:
                        idx1 = np.random.randint(0, self.swarmNum)
                        idx2 = np.random.randint(0, self.swarmNum)
                    if self.swarm[idx1][0] < self.swarm[idx2][0]:
                        newVel = self.PBest[idx1][2][port][2 * berth, :] - self.swarm[idx][2][port][2 * berth, :]
                    else:
                        newVel = self.PBest[idx2][2][port][2 * berth, :] - self.swarm[idx][2][port][2 * berth, :]
                newVel[newVel != 1] = 0
                newVel1 = c * np.random.rand(len(newVel)) * newVel
                newVel1[newVel1 > 1] = 1
                berthAllocationVel[port][2 * berth, :, 0] = newVel
                berthAllocationVel[port][2 * berth, :, 0] = newVel1
            for ves in range(len(self.ini.vesselsInPort[port])):
                Pc = np.random.rand()
                temp = 2 * berthNum + ves
                if Pc > self.Pc[idx]:
                    newVel = self.PBest[idx][2][port][temp, :] - self.swarm[idx][2][port][temp, :]
                else:
                    idx1 = np.random.randint(0, self.swarmNum)
                    idx2 = np.random.randint(0, self.swarmNum)
                    if idx1 == idx or idx2 == idx or idx1 == idx2:
                        idx1 = np.random.randint(0, self.swarmNum)
                        idx2 = np.random.randint(0, self.swarmNum)
                    if self.swarm[idx1][0] < self.swarm[idx2][0]:
                        newVel = self.PBest[idx1][2][port][temp, :] - self.swarm[idx][2][port][temp, :]
                    else:
                        newVel = self.PBest[idx2][2][port][temp, :] - self.swarm[idx][2][port][temp, :]
                newVel[newVel != 1] = 0
                newVel1 = c * np.random.rand(len(newVel)) * newVel
                newVel1[newVel1 > 1] = 1
                berthAllocationVel[port][temp, :, 0] = newVel
                berthAllocationVel[port][temp, :, 0] = newVel1
        return travelTimeVel, berthAllocationVel

    def positionSubtraction(self, c, position1, position2):
        travelTimeVel = []
        for ves in range(self.ini.vesselNum):
            travelTimeVel.append([])

        for ves in range(self.ini.vesselNum):
            for p in range(len(self.ini.portOrder[ves])):
                newVel = c * np.random.rand() * (position1[1][ves][p] - position2[1][ves][p])
                travelTimeVel[ves].append(newVel)

        berthAllocationVel = []
        for i in range(self.ini.portNum):
            temp = position1[2][i] - position2[2][i]
            temp[temp != 1] = 0
            temp1 = c * np.random.rand(temp.shape[0], temp.shape[1]) * temp
            temp1[temp1 > 1] = 1
            berthAllocationVel.append(np.stack((temp, temp1), 2))
        return travelTimeVel, berthAllocationVel

    def coefficientMultiplyVelocity(self, coefficient, velocity):
        travelTimeVelocity, berthAllocationVelocity = velocity

        for ves in range(self.ini.vesselNum):
            for p in range(len(self.ini.portOrder[ves])):
                travelTimeVelocity[ves][p] *= coefficient

        for i in range(self.ini.portNum):
            berthAllocationVelocity[i][:, :, 1] = berthAllocationVelocity[i][:, :, 1] * coefficient
        return travelTimeVelocity, berthAllocationVelocity

    def velocityAddition(self, velocity1, velocity2):
        travelTimeVel = []
        for ves in range(self.ini.vesselNum):
            travelTimeVel.append([])

        for ves in range(self.ini.vesselNum):
            for p in range(len(self.ini.portOrder[ves])):
                travelTimeVel[ves].append(velocity1[0][ves][p] + velocity2[0][ves][p])

        berthAllocationVel = []
        for i in range(self.ini.portNum):
            temp = velocity1[1][i][:, :, 0] + velocity2[1][i][:, :, 0]
            temp[temp >= 1] = 1
            berthAllocationVel.append(
                np.stack((temp, np.maximum(velocity1[1][i][:, :, 1], velocity2[1][i][:, :, 1])), 2))

        return travelTimeVel, berthAllocationVel

    def generateNewBerthAllocation(self, velocity, lastBerthAllocationArray, feasibleBerthAllocation):
        for p, port in enumerate(feasibleBerthAllocation):
            berthNum = self.ini.berthNum[p]
            vesselsInPort = self.ini.vesselsInPort[p]
            for b, berth in enumerate(port):
                if len(berth) > 0:
                    velocity[1][p][2 * b, :, :] = 0
                    lastBerthAllocationArray[p][2 * b, :] = 0
                    for v in range(len(berth) - 1):
                        idx = 2 * berthNum + vesselsInPort.index(berth[v])
                        nextIdx = 2 * berthNum + vesselsInPort.index(berth[v + 1])
                        velocity[1][p][:, idx, :] = 0
                        velocity[1][p][nextIdx, :, :] = 0
                        lastBerthAllocationArray[p][:, idx] = 0
                        lastBerthAllocationArray[p][nextIdx, :] = 0

        berthAllocationArray = []
        for order in range(self.ini.portNum):
            berthAllocationArray.append([])

        ports = np.arange(0, self.ini.portNum)
        np.random.shuffle(ports)

        for portCount, port in enumerate(ports):
            berthNum = self.ini.berthNum[port]
            vesselsInPort = self.ini.vesselsInPort[port]
            vesselNum = len(vesselsInPort)
            baShape = 2 * berthNum + vesselNum
            ba = np.zeros((baShape, baShape))

            alreadyArrangedVessels = []
            for b in feasibleBerthAllocation[port]:
                for v in b:
                    alreadyArrangedVessels.append(vesselsInPort.index(v))
            alreadyArrangedVessels = np.array(alreadyArrangedVessels)

            berths = np.arange(0, berthNum)
            np.random.shuffle(berths)

            for berth in berths:
                if len(feasibleBerthAllocation[port][berth]) == 0:
                    lastIdx = 2 * berth
                else:
                    lastIdx = feasibleBerthAllocation[port][berth][len(feasibleBerthAllocation[port][berth]) - 1]
                    lastIdx = 2 * berthNum + vesselsInPort.index(lastIdx)

                arr = velocity[1][port][lastIdx, 2 * berthNum:baShape, 1]
                arr = np.append(arr, velocity[1][port][lastIdx, 2 * berth + 1, 1])
                if np.max(arr) > 0:
                    probability = arr / np.sum(arr)
                    idx = np.random.choice(np.arange(0, vesselNum + 1), p=probability)
                    idx = 2 * berth + 1 if idx == vesselNum else idx + 2 * berthNum
                else:
                    arr = lastBerthAllocationArray[port][lastIdx, 2 * berthNum:baShape]
                    arr = np.append(arr, lastBerthAllocationArray[port][lastIdx, 2 * berth + 1])
                    if np.max(arr) > 0:
                        idx = np.argmax(arr)
                        idx = 2 * berth + 1 if idx == vesselNum else idx + 2 * berthNum
                    else:
                        finish = np.arange(vesselNum, vesselNum + berthNum)
                        arr = np.where(np.max(ba[:, berthNum * 2:baShape], axis=0) == 0)[0]
                        arr = np.setdiff1d(arr, alreadyArrangedVessels)
                        arr = np.append(arr, finish)
                        idx = np.random.choice(arr)
                        idx = 2 * berth + 1 if idx >= vesselNum else 2 * berthNum + idx
                ba[2 * berth, idx] = 1
                velocity[1][port][lastIdx, :, :] = 0
                velocity[1][port][:, idx, :] = 0
                lastBerthAllocationArray[port][lastIdx, :] = 0
                lastBerthAllocationArray[port][:, idx] = 0
                while idx != 2 * berth + 1:
                    arr = velocity[1][port][idx, 2 * berthNum:baShape, 1]
                    arr = np.append(arr, velocity[1][port][idx, 2 * berth + 1, 1])
                    if np.max(arr) != 0:
                        probability = arr / np.sum(arr)
                        nextIdx = np.random.choice(np.arange(0, vesselNum + 1), p=probability)
                        nextIdx = 2 * berth + 1 if nextIdx == vesselNum else nextIdx + 2 * berthNum
                    else:
                        arr = lastBerthAllocationArray[port][idx, 2 * berthNum:baShape]
                        arr = np.append(arr, lastBerthAllocationArray[port][idx, 2 * berth + 1])
                        if np.max(arr) > 0:
                            nextIdx = np.argmax(arr)
                            nextIdx = 2 * berth + 1 if nextIdx == vesselNum else nextIdx + 2 * berthNum
                        else:
                            finish = np.arange(vesselNum, vesselNum + berthNum)
                            arr = np.where(np.max(ba[:, 2 * berthNum:baShape], axis=0) == 0)[0]
                            arr = np.setdiff1d(arr, alreadyArrangedVessels)
                            arr = np.append(arr, finish)
                            nextIdx = np.random.choice(arr)
                            nextIdx = 2 * berth + 1 if nextIdx >= vesselNum else nextIdx + 2 * berthNum
                    ba[idx, nextIdx] = 1
                    velocity[1][port][idx, :, :] = 0
                    velocity[1][port][:, nextIdx, :] = 0
                    lastBerthAllocationArray[port][idx, :] = 0
                    lastBerthAllocationArray[port][:, nextIdx] = 0
                    idx = nextIdx
            arrangedVessels = np.where(np.max(ba[:, 2 * berthNum:baShape], axis=0) == 1)[0]
            arrangedVesselsIdx = arrangedVessels + 2 * berthNum
            totalArrangedVessels = np.append(arrangedVessels, alreadyArrangedVessels)
            notArrangedVessels = np.setdiff1d(np.arange(0, vesselNum), totalArrangedVessels)
            if len(notArrangedVessels) > 0:
                rList1 = []
                rList2 = []
                for row in range(0, 2 * berthNum, 2):
                    rList1.append(row)
                    rList2.append(row + 1)
                rList1 = np.array(rList1)
                rList2 = np.array(rList2)
                rList1 = np.concatenate((rList1, arrangedVesselsIdx))
                rList2 = np.concatenate((rList2, arrangedVesselsIdx))
                for vessel in notArrangedVessels:
                    if np.random.rand() > 0.5:
                        r = np.random.choice(rList1)
                        if np.max(ba[r, :]) > 0:
                            temp = np.argmax(ba[r, :])
                            ba[r, temp] = 0
                            ba[r, vessel + 2 * berthNum] = 1
                            ba[vessel + 2 * berthNum, temp] = 1
                        else:
                            ba[r, vessel + 2 * berthNum] = 1
                            ba[vessel + 2 * berthNum, r + 1] = 1
                        rList1 = np.append(rList1, vessel + 2 * berthNum)
                        rList2 = np.append(rList2, vessel + 2 * berthNum)
                    else:
                        r = np.random.choice(rList2)
                        if np.max(ba[:, r]) > 0:
                            temp = np.argmax(ba[:, r])
                            ba[temp, r] = 0
                            ba[vessel + 2 * berthNum, r] = 1
                            ba[temp, vessel + 2 * berthNum] = 1
                        else:
                            ba[vessel + 2 * berthNum, r] = 1
                            ba[r - 1, vessel + 2 * berthNum] = 1
                        rList1 = np.append(rList1, vessel + 2 * berthNum)
                        rList2 = np.append(rList2, vessel + 2 * berthNum)
            ba[0:2 * berthNum, 0:2 * berthNum] = 0
            berthAllocationArray[port] = ba
        berthAllocation = self.berthAllocationTransferToList(berthAllocationArray)
        return berthAllocation

    def positionUpdating(self):
        swarmVelocity = []
        for _ in range(self.swarmNum):
            swarmVelocity.append(self.initVelocity())
        flags = np.zeros(self.swarmNum)

        self.getPBest()
        self.getGBest()

        start = time.time()
        self.runtime = None
        sta = 0
        for iteration in range(self.iterationNum):
            omega = 0.9 - 0.5 * iteration / self.iterationNum
            swarm = []
            for solution in range(self.swarmNum):
                if flags[solution] < self.rg:
                    velocity1 = self.coefficientMultiplyVelocity(omega, swarmVelocity[solution])
                    velocity2 = self.CLPSOPositionSubtraction(self.c1, solution)
                    velocity = self.velocityAddition(velocity1, velocity2)
                    for ves in range(self.ini.vesselNum):
                        for p in range(len(self.ini.portOrder[ves])):
                            if velocity[0][ves][p] < -self.travelTimeVelocityLimit[ves][p]:
                                velocity[0][ves][p] = -self.travelTimeVelocityLimit[ves][p]
                            elif velocity[0][ves][p] > self.travelTimeVelocityLimit[ves][p]:
                                velocity[0][ves][p] = self.travelTimeVelocityLimit[ves][p]
                    swarmVelocity[solution] = copy.deepcopy(velocity)
                else:
                    velocity1 = self.coefficientMultiplyVelocity(omega, swarmVelocity[solution])
                    velocity2 = self.positionSubtraction(self.c1, self.PBest[solution], self.swarm[solution])
                    velocity3 = self.positionSubtraction(self.c2, self.GBest, self.swarm[solution])
                    velocity = self.velocityAddition(velocity1, velocity2)
                    velocity = self.velocityAddition(velocity, velocity3)
                    for ves in range(self.ini.vesselNum):
                        for p in range(len(self.ini.portOrder[ves])):
                            if velocity[0][ves][p] < -self.travelTimeVelocityLimit[ves][p]:
                                velocity[0][ves][p] = -self.travelTimeVelocityLimit[ves][p]
                            elif velocity[0][ves][p] > self.travelTimeVelocityLimit[ves][p]:
                                velocity[0][ves][p] = self.travelTimeVelocityLimit[ves][p]
                    swarmVelocity[solution] = copy.deepcopy(velocity)
                    flags[solution] = 0

                travelTimeVel = velocity[0]
                newTravelTime = []

                for ves in range(self.ini.vesselNum):
                    newTravelTime.append([])

                for ves in range(self.ini.vesselNum):
                    for p in range(len(self.ini.portOrder[ves])):
                        t = self.swarm[solution][1][ves][p] + travelTimeVel[ves][p]
                        if t > self.ini.travelTimeLimit[ves][p, 1] or t < self.ini.travelTimeLimit[ves][p, 0]:
                            newTravelTime = copy.deepcopy(self.swarm[solution][1])
                            break
                        else:
                            newTravelTime[ves].append(t)
                    else:
                        continue
                    break

                alpha = np.random.rand()

                for p in range(self.ini.portNum):
                    velocity[1][p][:, :, 0][velocity[1][p][:, :, 1] < alpha] = 0
                    velocity[1][p][:, :, 1][velocity[1][p][:, :, 1] < alpha] = 0

                portOrder = copy.deepcopy(self.ini.portOrder)
                berthAllocation = []
                for order in range(self.ini.portNum):
                    ba = []
                    for j in range(self.ini.berthNum[order]):
                        ba.append([])
                    berthAllocation.append(ba)

                tryBerthAllocation = self.generateNewBerthAllocation(velocity, self.swarm[solution][2], berthAllocation)

                sumLen1 = 0
                for order in portOrder:
                    sumLen1 += len(order)
                while sumLen1 != 0:
                    for vessel, order in enumerate(portOrder):
                        if order.shape[0] > 0:
                            for berth, berthOrder in enumerate(tryBerthAllocation[order[0]]):
                                if len(berthOrder) != 0 and vessel == berthOrder[0]:
                                    portOrder[vessel] = np.delete(portOrder[vessel], 0)
                                    berthAllocation[order[0]][berth].append(tryBerthAllocation[order[0]][berth][0])
                                    del tryBerthAllocation[order[0]][berth][0]
                                    break
                    sumLen2 = 0
                    for order in portOrder:
                        sumLen2 += len(order)
                    if sumLen1 == sumLen2:
                        tryBerthAllocation = self.generateNewBerthAllocation(velocity, self.swarm[solution][2],
                                                                             berthAllocation)
                    sumLen1 = sumLen2
                totalCost = self.ini.calculateObjectiveFunction(newTravelTime, berthAllocation)
                if totalCost < self.PBest[solution][0]:
                    flags[solution] = 0
                else:
                    flags[solution] += 1
                swarm.append([totalCost, newTravelTime, self.berthAllocationTransferToArray(berthAllocation)])
            self.swarm = swarm
            previousCost = self.GBest[0]
            self.getPBest()
            self.getGBest()
            print(str(iteration) + ' ' + str(self.GBest[0]))
            print(self.GBest[1])
            print(self.berthAllocationTransferToList(self.GBest[2]))
            if (previousCost - self.GBest[0]) / previousCost < 0.00001:
                sta += 1
            else:
                sta = 0
            print('sta:' + str(sta))
            if sta == 100:
                self.runtime = time.time() - start
            if iteration == self.iterationNum - 1:
                if self.runtime is None:
                    self.runtime = time.time() - start
        with open('PSOresult.txt', 'a') as f:
            f.write('PSO_portNum{}_berthNum{}_vesselNum{}\n'.format(self.ini.portNum, self.ini.berthNum, self.ini.vesselNum))
            f.write('totalCost:{}\n'.format(self.GBest[0]))
            f.write('berthAllocation:{}\n'.format(self.berthAllocationTransferToList(self.GBest[2])))
            travelTime = str(self.GBest[1]).replace('array', 'np.array')
            f.write('travelTime:{}\n'.format(travelTime))
            f.write('runTime:{}\n'.format(self.runtime))

    def run(self):
        try:
            self.positionUpdating()
        except BaseException as e:
            if isinstance(e, KeyboardInterrupt):
                with open('PSOresult.txt', 'a') as f:
                    f.write('PSO_portNum{}_berthNum{}_vesselNum{}\n'.format(self.ini.portNum, self.ini.berthNum, self.ini.vesselNum))
                    f.write('totalCost:{}\n'.format(self.GBest[0]))
                    f.write('berthAllocation:{}\n'.format(self.berthAllocationTransferToList(self.GBest[2])))
                    travelTime = str(self.GBest[1]).replace('array', 'np.array')
                    f.write('travelTime:{}\n'.format(travelTime))
                    f.write('runTime:{}\n'.format(self.runtime))


if __name__ == '__main__':
    ini = Initialization(portNum=4, berthNum=4, vesselNum=5)
    pso = PSO(ini)
    pso.run()
    # pso.ini.paintSolution(pso.GBest[1], pso.berthAllocationTransferToList(pso.GBest[2]))
