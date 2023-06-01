import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib


class Initialization:
    def __init__(self, portNum, berthNum, vesselNum):
        np.random.seed(2)
        self.portNum = portNum
        self.berthNum = np.random.randint(berthNum, berthNum + 1, self.portNum)
        self.vesselNum = vesselNum
        self.handlingTimes = self.initHandlingTime()
        self.distances = self.initDistance()
        print(self.distances)
        self.portOrder, self.vesselsInPort = self.initPortOrder()
        self.speedLimit = [14, 19]
        self.travelTimeLimit = self.initTravelTimeLimit()
        self.timeWindows = self.initTimeWindows()
        np.random.seed(None)

        with open('ini.txt', 'w') as f:
            f.write('portNum:{}\n'.format(portNum))
            f.write('berthNum:\n')
            for b in self.berthNum:
                f.write(str(b) + ',')
            f.write('\n')
            f.write('vesselNum:{}\n'.format(vesselNum))
            f.write('handlingTimes:\n')
            for p in range(self.portNum):
                for b in range(self.berthNum[p]):
                    for v in range(self.vesselNum):
                        f.write(str(self.handlingTimes[p][b, v]) + ',')
                    f.write('\n')
                f.write('\n')
            f.write('distances1:\n')
            for p1 in range(self.portNum):
                for p2 in range(self.portNum):
                    f.write(str(self.distances[0][p1, p2]) + ',')
                f.write('\n')
            f.write('distances2:\n')
            for v in range(self.vesselNum):
                f.write(str(self.distances[1][v]) + ',')
            f.write('\n')
            f.write('portOrder:\n')
            for v in range(self.vesselNum):
                for p in self.portOrder[v]:
                    f.write(str(self.portOrder[v][p]) + ',')
                f.write('\n')
            f.write('timeWindows:\n')
            for v in range(self.vesselNum):
                for i in range(2):
                    for p in range(len(self.portOrder[v])):
                        f.write(str(self.timeWindows[v][p, i]) + ',')
                    f.write('\n')
                f.write('\n')

        self.Fc = 250
        self.fd = 42 / 24
        self.Hc = 200
        self.Ic = 200
        self.Dc = 300

    def initHandlingTime(self):
        handlingTimes = []
        for i in range(self.portNum):
            handlingTimes.append(np.random.uniform(3, 15, (self.berthNum[i], self.vesselNum)))
        return handlingTimes

    def initDistance(self):
        portDistances = np.random.uniform(160, 200, (self.portNum, self.portNum))
        portDistances = np.triu(portDistances)
        portDistances += portDistances.T - np.diag(portDistances.diagonal())
        vesselDistances = np.random.uniform(160, 200, self.vesselNum)
        return portDistances, vesselDistances

    def initPortOrder(self):
        portOrder = []
        for i in range(self.vesselNum):
            arr = np.arange(0, self.portNum)
            # arr = np.random.choice(arr, np.random.randint(2, self.portNum + 1), replace=False)
            portOrder.append(arr)

        vesselsInPort = []
        for i in range(self.portNum):
            vesselsInPort.append([])
        for i in range(len(portOrder)):
            for j in portOrder[i]:
                vesselsInPort[j].append(i)
        return portOrder, vesselsInPort

    def initTravelTimeLimit(self):
        travelTimeLimit = []
        for v in range(self.vesselNum):
            arr = np.zeros((len(self.portOrder[v]), 2))
            for p, port in enumerate(self.portOrder[v]):
                if p == 0:
                    arr[p, 0] = self.distances[1][v] / self.speedLimit[1]
                    arr[p, 1] = self.distances[1][v] / self.speedLimit[0]
                if p != len(self.portOrder[v]) - 1:
                    arr[p + 1, 0] = self.distances[0][port, self.portOrder[v][p + 1]] / self.speedLimit[1]
                    arr[p + 1, 1] = self.distances[0][port, self.portOrder[v][p + 1]] / self.speedLimit[0]
            travelTimeLimit.append(arr)
        return travelTimeLimit

    def initTimeWindows(self):
        timeWindows = []
        for v in range(self.vesselNum):
            arr = np.zeros((len(self.portOrder[v]), 2))
            for p, port in enumerate(self.portOrder[v]):
                if p == 0:
                    arr[p, 0] = np.random.uniform(self.travelTimeLimit[v][p, 0], self.travelTimeLimit[v][p, 1])
                    arr[p, 1] = arr[p, 0] + np.random.uniform(3, 15)
                if p != len(self.portOrder[v]) - 1:
                    arr[p + 1, 0] = arr[p, 1] + np.random.uniform(self.travelTimeLimit[v][p + 1, 0],
                                                                  self.travelTimeLimit[v][p + 1, 1])
                    arr[p + 1, 1] = arr[p + 1, 0] + np.random.uniform(3, 15)
            timeWindows.append(arr)
        return timeWindows

    def generateFeasibleSolution(self):
        travelTime = []
        for v in range(self.vesselNum):
            temp = []
            for p in range(len(self.portOrder[v])):
                temp.append(np.random.uniform(self.travelTimeLimit[v][p, 0], self.travelTimeLimit[v][p, 1]))
            travelTime.append(np.array(temp))

        tryBerthAllocation = []
        for order in range(self.portNum):
            arr = np.random.choice(self.vesselsInPort[order], len(self.vesselsInPort[order]), replace=False)
            ba = []
            for j in range(self.berthNum[order]):
                ba.append([])
            for j in arr:
                r = np.random.randint(0, self.berthNum[order])
                ba[r].append(j)
            tryBerthAllocation.append(ba)

        totalCost, berthAllocation = self.calculateObjectiveFunction(travelTime, tryBerthAllocation,
                                                                     isGenerateFeasibleSolution=True)

        return totalCost, travelTime, berthAllocation

    def calculateObjectiveFunction(self, travelTime, tryBerthAllocation, isGenerateFeasibleSolution=False,
                                   outputTime=False, handlingBeginTimeTable=None):
        flag = False
        if handlingBeginTimeTable is None:
            flag = True

        arriveTable = []
        handlingEndTimeTable = []
        deltaEFTTable = []
        for order in range(self.vesselNum):
            arriveTable.append([])
            handlingEndTimeTable.append([])
            deltaEFTTable.append([])

        if flag:
            handlingBeginTimeTable = []
            for order in range(self.vesselNum):
                handlingBeginTimeTable.append([])

        berthTimeTable = []
        for order in range(self.portNum):
            temp = []
            for j in range(self.berthNum[order]):
                temp.append(0)
            berthTimeTable.append(temp)

        tryBerthAllocation = copy.deepcopy(tryBerthAllocation)
        berthAllocation = None
        if isGenerateFeasibleSolution:
            berthAllocation = []
            for order in range(self.portNum):
                ba = []
                for j in range(self.berthNum[order]):
                    ba.append([])
                berthAllocation.append(ba)

        travelTime1 = copy.deepcopy(travelTime)
        travelTimeToFirstPort = []
        for order, time in enumerate(travelTime1):
            travelTimeToFirstPort.append(time[0])
            travelTime1[order] = np.delete(travelTime1[order], 0)
        travelTimeToFirstPort = np.array(travelTimeToFirstPort)
        tempArriveTable = travelTimeToFirstPort
        for order in range(self.vesselNum):
            arriveTable[order].append(tempArriveTable[order])

        portOrder = copy.deepcopy(self.portOrder)
        timeWindows = copy.deepcopy(self.timeWindows)

        sumLen1 = 0
        for order in portOrder:
            sumLen1 += len(order)
        while sumLen1 != 0:
            for vessel, order in enumerate(portOrder):
                if order.shape[0] > 0:
                    for berth, berthOrder in enumerate(tryBerthAllocation[order[0]]):
                        if len(berthOrder) != 0 and vessel == berthOrder[0]:
                            if flag:
                                handlingBeginTime = max(timeWindows[vessel][0, 0], tempArriveTable[vessel],
                                                        berthTimeTable[order[0]][berth])
                                handlingBeginTimeTable[vessel].append(handlingBeginTime)

                            handlingBeginTime = handlingBeginTimeTable[vessel][
                                np.where(self.portOrder[vessel] == order[0])[0][0]]
                            handlingTime = self.handlingTimes[order[0]][berth, vessel]
                            handlingEndTime = handlingBeginTime + handlingTime
                            handlingEndTimeTable[vessel].append(handlingEndTime)

                            if handlingEndTime - timeWindows[vessel][0, 1] > 0:
                                deltaEFTTable[vessel].append(handlingEndTime - timeWindows[vessel][0, 1])
                            else:
                                deltaEFTTable[vessel].append(0)

                            berthTimeTable[order[0]][berth] = handlingEndTime
                            if order.shape[0] > 1:
                                nAT = handlingEndTime + travelTime1[vessel][0]
                                tempArriveTable[vessel] = nAT
                                arriveTable[vessel].append(nAT)

                            portOrder[vessel] = np.delete(portOrder[vessel], 0)
                            if isGenerateFeasibleSolution:
                                berthAllocation[order[0]][berth].append(tryBerthAllocation[order[0]][berth][0])
                            del tryBerthAllocation[order[0]][berth][0]
                            timeWindows[vessel] = np.delete(timeWindows[vessel], [0], 0)
                            if travelTime1[vessel].shape[0] != 0:
                                travelTime1[vessel] = np.delete(travelTime1[vessel], 0)
                            break

            sumLen2 = 0
            for order in portOrder:
                sumLen2 += len(order)
            if sumLen2 == sumLen1:
                if isGenerateFeasibleSolution:
                    for i, port in enumerate(tryBerthAllocation):
                        vessel = []
                        for berth in port:
                            vessel.extend(berth)
                        if len(vessel) > 0:
                            arr = np.random.choice(vessel, len(vessel), replace=False)
                            tryBerthAllocation[i] = []
                            for j in range(self.berthNum[i]):
                                tryBerthAllocation[i].append([])
                            for j in arr:
                                r = np.random.randint(0, self.berthNum[i])
                                tryBerthAllocation[i][r].append(j)
                else:
                    return None
            sumLen1 = sumLen2

        arriveList = []
        handlingBeginTimeList = []
        handlingEndTimeList = []
        deltaEFTList = []
        for i in range(self.vesselNum):
            arriveList.extend(arriveTable[i])
            handlingBeginTimeList.extend(handlingBeginTimeTable[i])
            handlingEndTimeList.extend(handlingEndTimeTable[i])
            deltaEFTList.extend(deltaEFTTable[i])
        arriveList = np.array(arriveList)
        handlingBeginTimeList = np.array(handlingBeginTimeList)
        handlingEndTimeList = np.array(handlingEndTimeList)
        deltaEFTList = np.array(deltaEFTList)

        idleCost = np.sum(handlingBeginTimeList - arriveList) * self.Ic
        handlingCost = np.sum(handlingEndTimeList - handlingBeginTimeList) * self.Hc
        delayCost = np.sum(deltaEFTList) * self.Dc

        fuelConsumptionCost = 0
        for v, tt in enumerate(travelTime):
            for p, t in enumerate(tt):
                if p == 0:
                    temp = self.Fc * self.fd * np.power(self.distances[1][v] / self.speedLimit[0], 3) * np.power(t, -2)
                else:
                    lastPort = self.portOrder[v][p - 1]
                    port = self.portOrder[v][p]
                    temp = self.Fc * self.fd * np.power(self.distances[0][lastPort, port] / self.speedLimit[0],
                                                        3) * np.power(t, -2)
                fuelConsumptionCost += temp
        totalCost = idleCost + handlingCost + delayCost + fuelConsumptionCost

        print('idleCost:' + str(idleCost))
        print('handlingCost' + str(handlingCost))
        print('delayCost' + str(delayCost))
        print('fuelConsumptionCost' + str(fuelConsumptionCost))
        print('totalCost' + str(totalCost))

        if isGenerateFeasibleSolution:
            return totalCost, berthAllocation
        else:
            if outputTime:
                return totalCost, travelTime, arriveTable, handlingBeginTimeTable, handlingEndTimeTable
            return totalCost

    def paintSolution(self, travelTime, berthAllocation, handlingBeginTimeTable=None):
        totalCost, travelTime, arriveTable, handlingBeginTimeTable, handlingEndTimeTable = self.calculateObjectiveFunction(
            travelTime, berthAllocation, outputTime=True, handlingBeginTimeTable=handlingBeginTimeTable)

        print(arriveTable)
        print(handlingBeginTimeTable)
        print(handlingEndTimeTable)
        cnames = ['black','black','black','black','black','black', 'dimgray', 'gray', 'darkgray', 'silver', 'lightgray', 'darkviolet', 'chocolate', 'darkblue', 'tan',
                  'deepskyblue', 'orange', 'purple', 'red', 'coral', 'cornflowerblue', 'darkcyan', 'darkgoldenrod',
                  'olive', 'darkgreen', 'darksalmon', 'darkseagreen', 'darkslateblue', 'darkslategray',
                  'darkturquoise', 'green', 'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange',
                  'darkorchid', 'cyan', 'deeppink', 'blue', 'dodgerblue', 'firebrick', 'floralwhite', 'forestgreen',
                  'fuchsia', 'gainsboro', 'gold', 'goldenrod', 'violet', 'greenyellow', 'honeydew', 'hotpink',
                  'indianred', 'indigo', 'ivory', 'khaki', 'lavender', 'lavenderblush', 'lawngreen', 'lemonchiffon',
                  'lime', 'limegreen', 'linen', 'magenta', 'maroon', 'mediumaquamarine', 'mediumblue', 'mediumorchid',
                  'mediumpurple', 'mediumseagreen', 'mediumslateblue', 'mediumspringgreen', 'mediumturquoise',
                  'mediumvioletred', 'midnightblue', 'mintcream', 'mistyrose', 'moccasin', 'navajowhite', 'navy',
                  'oldlace', 'olivedrab', 'orangered', 'orchid', 'palegoldenrod', 'palegreen', 'paleturquoise',
                  'palevioletred', 'papayawhip', 'peachpuff', 'peru', 'pink', 'plum', 'powderblue', 'rosybrown',
                  'royalblue', 'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 'seashell', 'sienna', 'brown',
                  'skyblue', 'slateblue', 'slategray', 'snow', 'springgreen', 'steelblue', 'teal', 'thistle', 'tomato',
                  'turquoise', 'wheat', 'white', 'whitesmoke', 'yellow', 'yellowgreen', 'aliceblue', 'antiquewhite',
                  'azure', 'beige', 'blanchedalmond', 'bisque', 'cornsilk', 'aquamarine', 'blueviolet', 'cadetblue',
                  'aqua', 'ghostwhite', 'lightblue', 'lightcoral', 'lightcyan', 'lightgoldenrodyellow', 'lightgreen',
                  'darkred', 'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightslategray',
                  'lightsteelblue', 'lightyellow', 'chartreuse', 'burlywood']

        del matplotlib.font_manager.weight_dict['roman']
        matplotlib.font_manager._rebuild()
        plt.rc('font', family='Times New Roman')
        for port in range(self.portNum):
            firstPort = []
            speedsToNextPort = []
            nextPort = []
            arriveTimeInPort = []
            handlingBeginTimeInPort = []
            handlingEndTimeInPort = []
            timeWindowInPort = []
            for v, order in enumerate(self.portOrder):
                if port == order[0]:
                    firstPort.append([np.round(self.distances[1][v] / travelTime[v][0], decimals=1)])
                else:
                    firstPort.append([])
                idx = np.where(order == port)[0]
                if len(idx) > 0:
                    idx = idx[0]
                    if idx + 1 != len(order):
                        speedsToNextPort.append([np.round(
                            self.distances[0][port, order[idx + 1]] / travelTime[v][idx + 1], decimals=1)])
                        nextPort.append([order[idx + 1]])
                    else:
                        speedsToNextPort.append([])
                        nextPort.append([])
                    arriveTimeInPort.append([arriveTable[v][idx]])
                    handlingBeginTimeInPort.append([handlingBeginTimeTable[v][idx]])
                    handlingEndTimeInPort.append([handlingEndTimeTable[v][idx]])
                    timeWindowInPort.append(self.timeWindows[v][idx])
                else:
                    speedsToNextPort.append([])
                    nextPort.append([])
                    arriveTimeInPort.append([])
                    handlingBeginTimeInPort.append([])
                    handlingEndTimeInPort.append([])
                    timeWindowInPort.append([])

            maxHandlingEndTimeInPort = 0
            minArriveTimeInPort = np.inf
            maxEFTInPort = 0
            minStartInPort = np.inf
            for v in range(self.vesselNum):
                if len(handlingEndTimeInPort[v]) > 0:
                    if maxHandlingEndTimeInPort < handlingEndTimeInPort[v][0]:
                        maxHandlingEndTimeInPort = handlingEndTimeInPort[v][0]
                if len(arriveTimeInPort[v]) > 0:
                    if minArriveTimeInPort > arriveTimeInPort[v][0]:
                        minArriveTimeInPort = arriveTimeInPort[v][0]
                if len(timeWindowInPort[v]) > 0:
                    if timeWindowInPort[v][0] < minStartInPort:
                        minStartInPort = timeWindowInPort[v][0]
                    if timeWindowInPort[v][1] > maxEFTInPort:
                        maxEFTInPort = timeWindowInPort[v][1]
            upperLimit = max(maxHandlingEndTimeInPort, maxEFTInPort)
            lowerLimit = min(minStartInPort, minArriveTimeInPort)
            upperLimit = int(np.ceil(upperLimit)) + 1
            lowerLimit = int(np.floor(lowerLimit)) - 1

            col = (upperLimit - lowerLimit) / 5
            row = self.berthNum[port] * 3
            plt.figure(dpi=200, figsize=(row, col))
            ax = plt.axes([0.1, 0.1, 0.8, 0.8])
            plt.xticks([x for x in range(self.berthNum[port] + 1)], color='w')
            plt.yticks([y for y in range(lowerLimit, upperLimit, 1)])
            plt.ylim(lowerLimit, upperLimit)
            ax.set_title('port' + str(port+1))

            for b in range(self.berthNum[port]):
                plt.text(0.1 + 1 / self.berthNum[port] * b, -0.03, 'berth' + str(b+1), transform=ax.transAxes)

            temp1 = 1 / self.berthNum[port]
            temp2 = temp1 * 0.2
            for v in range(self.vesselNum):
                if len(arriveTimeInPort[v]) > 0:
                    bForV = None
                    for b, berth in enumerate(berthAllocation[port]):
                        if v in berth:
                            bForV = b
                            break
                    plt.axhline(y=arriveTimeInPort[v][0], xmin=bForV * temp1 + temp2, xmax=(bForV + 1) * temp1 - temp2,
                                c=cnames[v], ls="-", lw=0.5, marker='o', markersize=3)
                    plt.axhline(y=timeWindowInPort[v][0], xmin=bForV * temp1 + temp2, xmax=(bForV + 1) * temp1 - temp2,
                                c=cnames[v], ls="-.", lw=0.5, marker='v', markersize=3)
                    plt.axhline(y=timeWindowInPort[v][1], xmin=bForV * temp1 + temp2, xmax=(bForV + 1) * temp1 - temp2,
                                c=cnames[v], ls="--", lw=0.5, marker='s', markersize=3)
                    heightOfRect = handlingEndTimeInPort[v][0] - handlingBeginTimeInPort[v][0]
                    rect = plt.Rectangle(xy=(bForV + 0.2, handlingBeginTimeInPort[v][0]), width=0.6,
                                         height=heightOfRect, fill=False, ec=cnames[v])
                    ax.add_patch(rect)
                    plt.text(bForV + 0.3, handlingBeginTimeInPort[v][0] + heightOfRect / 1.8, 'vessel' + str(v+1),
                             fontsize=10, color=cnames[v])

                    if len(firstPort[v]) != 0:
                        plt.text(bForV + 0.05, handlingBeginTimeInPort[v][0] + heightOfRect / 2.2, str(firstPort[v][0]),
                                 fontsize=10, color=cnames[v])
                        ax.arrow(bForV + 0.125, handlingBeginTimeInPort[v][0] + heightOfRect / 2.2 - 1, 0, 1,
                                 length_includes_head=True, head_length=0.5, head_width=0.05, fc=cnames[v],
                                 ec=cnames[v])

                    if len(speedsToNextPort[v]) != 0:
                        plt.text(bForV + 0.85, handlingBeginTimeInPort[v][0] + heightOfRect / 2.2,
                                 str(speedsToNextPort[v][0]), fontsize=10, color=cnames[v])
                        ax.arrow(bForV + 0.925, handlingBeginTimeInPort[v][0] + heightOfRect / 2.2 + 0.8, 0, 1,
                                 length_includes_head=True, head_length=0.5, head_width=0.05, fc=cnames[v],
                                 ec=cnames[v])
                        # plt.text(bForV + 0.85, handlingBeginTimeInPort[v][0] + heightOfRect / 2.2 + 1.8,
                        #          'port' + str(nextPort[v][0]), fontsize=10, color=cnames[v])

            plt.savefig('port' + str(port) + '.pdf')


if __name__ == '__main__':
    ini = Initialization(portNum=3, berthNum=3, vesselNum=6)
    ini.paintSolution([[4.769215, 4.363433, 7.0557346],
                       [5.1229234, 4.363293, 7.0588],
                       [5.782788, 4.363427, 7.058687],
                       [5.3044424, 4.363425, 7.0586863],
                       [4.0182056, 3.8135436, 6.431334],
                       [7.0947237, 4.363435, 7.0588]],
                      [[[4, 1], [2, 3], [0, 5]],
                       [[4, 3], [0, 2], [1, 5]],
                       [[2, 5], [4, 1], [0, 3]]])
    # t, s, b = ini.generateFeasibleSolution()
    # print(t)
    # print(s)
    # print(b)
    # print(ini.timeWindows)
    # print(ini.handlingTimes)

    # ini.calculateObjectiveFunction([[12.627814, 13.318518, 13.911149, 14.052738],
    #                                 [12.971975, 11.977255, 13.509168, 14.052738],
    #                                 [13.166245, 13.3181715, 13.91117, 14.051858],
    #                                 [13.789268, 13.316363, 13.910962, 14.052375],
    #                                 [10.760545, 11.887967, 13.910993, 14.052734],
    #                                 [10.920155, 13.318482, 13.910922, 14.052738]],
    #                                [[[4], [3], [0, 2], [5, 1]],
    #                                 [[4, 0], [1], [5, 2], [3]],
    #                                 [[0], [1, 3, 2], [4], [5]],
    #                                 [[5, 0], [3], [2], [1, 4]]])
