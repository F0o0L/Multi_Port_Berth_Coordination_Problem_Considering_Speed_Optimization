package CLPSO_MVNS;

import CLPSO_MVNS.MyLinkedList.Node;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
//import java.util.LinkedList;
import java.util.Random;
import java.util.TreeMap;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;

public class CLPSO_MVNS {
	int swarmNum;
	int iterationNum;
	int PSOiterationNum;
	int VNSiterationNum;
	int portNum;
	int[] berthNum;
	int vesselNum;
	float[][][] handlingTimes;
	float[][] distances1;
	float[] distances2;
	int[][] portOrder;
	ArrayList<ArrayList<Integer>> vesselsInPort;
	float[][][] timeWindows;
	float[][][] travelTimeLimit;
	int[] speedLimit = { 14, 19 };
	float Fc = 250;
	float fd = 42.0f / 24;
	int Hc = 200;
	int Ic = 200;
	int Dc = 500;

	float[][] travelTimeVolocityLimit;
	int rg = 7;
	float[] pc;
	float c1 = 2;
	float c2 = 2;

	float[][][] travelTimeSwarm;
	float[] travelTimeSwarmCost;
	float[][][] travelTimeSwarmVelocity;
	ArrayList<ArrayList<MyLinkedList<Integer>>> berthAllocation;
	float[] PBestCost;
	float[][][] PBest;
	float GBestCost;
	float[][] GBest;
	int[] flags;
	
	float VNSLastCost;
	float VNSBestCost;
	ArrayList<ArrayList<MyLinkedList<Integer>>> VNSBestBerthAllocation;
	
	float bestCost;
	float[][] bestTravelTime;
	ArrayList<ArrayList<MyLinkedList<Integer>>> bestBerthAllocation;

	public CLPSO_MVNS(int swarmNum, int iterationNum, String fileName) {
		this.swarmNum = swarmNum;
		this.iterationNum = iterationNum;
		readFile(fileName);
	}
	
	public CLPSO_MVNS(int swarmNum, int iterationNum, int portNum,int berthNum,int vesselNum) {
		this.swarmNum = swarmNum;
		this.iterationNum = iterationNum;
		generateParameters(portNum, berthNum, vesselNum);
	}

	public void generateParameters(int portNum, int berthNum, int vesselNum) {
		Random random = new Random(5);
		this.portNum = portNum;
		this.berthNum = random.ints(portNum, berthNum, berthNum + 1).toArray();
		this.vesselNum = vesselNum;

		RandomGenerator rg = RandomGeneratorFactory.createRandomGenerator(random);
		rg.setSeed(5);
		UniformRealDistribution urd0 = new UniformRealDistribution(rg, 3, 15);
		handlingTimes = new float[portNum][][];
		for (int i = 0; i < portNum; i++) {
			handlingTimes[i] = new float[this.berthNum[i]][vesselNum];
			for (int j = 0; j < this.berthNum[i]; j++) {
				for (int k = 0; k < vesselNum; k++) {
					handlingTimes[i][j][k] = (float) urd0.sample();
				}
			}
		}
		UniformRealDistribution urd1 = new UniformRealDistribution(rg, 60, 100);
		distances1 = new float[portNum][portNum];
		for (int i = 0; i < portNum; i++) {
			for (int j = 0; j < portNum; j++) {
				distances1[i][j] = (float) urd1.sample();
			}
		}
		distances2 = new float[vesselNum];
		for (int i = 0; i < vesselNum; i++) {
			distances2[i] = (float) urd1.sample();
		}
		UniformIntegerDistribution uid = new UniformIntegerDistribution(rg, portNum, portNum);
		int[] temp = new int[portNum];
		for (int i = 0; i < portNum; i++) {
			temp[i] = i;
		}
		this.portOrder = new int[vesselNum][];
		for (int i = 0; i < vesselNum; i++) {
//			ArrayUtils.shuffle(temp, random);
			int pl = uid.sample();
			portOrder[i] = new int[pl];
			for (int j = 0; j < portOrder[i].length; j++) {
				portOrder[i][j] = temp[j];
			}
		}
		getExtraParameters();
		timeWindows = new float[vesselNum][][];
		for (int i = 0; i < vesselNum; i++) {
			timeWindows[i] = new float[portOrder[i].length][2];
			for (int j = 0; j < portOrder[i].length - 1; j++) {
				if (j == 0) {
					timeWindows[i][j][0] = (float) new UniformRealDistribution(rg, travelTimeLimit[i][0][0],
							travelTimeLimit[i][0][1]).sample();
					timeWindows[i][j][1] = timeWindows[i][j][0] + (float) urd0.sample();
				}
				timeWindows[i][j + 1][0] = timeWindows[i][j][1] + (float) new UniformRealDistribution(rg,
						travelTimeLimit[i][j + 1][0], travelTimeLimit[i][j + 1][1]).sample();
				timeWindows[i][j + 1][1] = timeWindows[i][j + 1][0] + (float) urd0.sample();
			}
		}
	}

	public void readFile(String fileName) {
		try {
			FileReader fr = new FileReader(fileName);
			BufferedReader br = new BufferedReader(fr);
			String[] temp = br.readLine().split(":");
			portNum = Integer.parseInt(temp[1]);
			br.readLine();
			temp = br.readLine().split(",");
			berthNum = new int[portNum];
			for (int i = 0; i < portNum; i++) {
				berthNum[i] = Integer.parseInt(temp[i]);
			}
			temp = br.readLine().split(":");
			vesselNum = Integer.parseInt(temp[1]);
			br.readLine();
			handlingTimes = new float[portNum][][];
			for (int i = 0; i < portNum; i++) {
				handlingTimes[i] = new float[berthNum[i]][vesselNum];
				for (int j = 0; j < berthNum[i]; j++) {
					temp = br.readLine().split(",");
					for (int k = 0; k < vesselNum; k++) {
						handlingTimes[i][j][k] = Float.parseFloat(temp[k]);
					}
				}
				br.readLine();
			}
			br.readLine();
			distances1 = new float[portNum][portNum];
			for (int i = 0; i < portNum; i++) {
				temp = br.readLine().split(",");
				for (int j = 0; j < portNum; j++) {
					distances1[i][j] = Float.parseFloat(temp[j]);
				}
			}
			br.readLine();
			distances2 = new float[vesselNum];
			temp = br.readLine().split(",");
			for (int i = 0; i < vesselNum; i++) {
				distances2[i] = Float.parseFloat(temp[i]);
			}
			br.readLine();
			this.portOrder = new int[vesselNum][];
			for (int i = 0; i < vesselNum; i++) {
				temp = br.readLine().split(",");
				portOrder[i] = new int[temp.length];
				for (int j = 0; j < temp.length; j++) {
					portOrder[i][j] = Integer.parseInt(temp[j]);
				}
			}
			br.readLine();
			timeWindows = new float[vesselNum][][];
			for (int i = 0; i < vesselNum; i++) {
				timeWindows[i] = new float[portOrder[i].length][2];
				for (int j = 0; j < 2; j++) {
					temp = br.readLine().split(",");
					for (int k = 0; k < portOrder[i].length; k++) {
						timeWindows[i][k][j] = Float.parseFloat(temp[k]);
					}
				}
				br.readLine();
			}
			br.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		getExtraParameters();
	}

	public void getExtraParameters() {
		travelTimeLimit = new float[vesselNum][][];
		for (int i = 0; i < vesselNum; i++) {
			travelTimeLimit[i] = new float[portOrder[i].length][2];
			for (int j = 0; j < portOrder[i].length; j++) {
				if (j == 0) {
					travelTimeLimit[i][j][0] = distances2[i] / speedLimit[1];
					travelTimeLimit[i][j][1] = distances2[i] / speedLimit[0];
				} else {
					travelTimeLimit[i][j][0] = distances1[portOrder[i][j - 1]][portOrder[i][j]] / speedLimit[1];
					travelTimeLimit[i][j][1] = distances1[portOrder[i][j - 1]][portOrder[i][j]] / speedLimit[0];
				}
			}
		}
		vesselsInPort = new ArrayList<ArrayList<Integer>>();
		for (int i = 0; i < portNum; i++) {
			vesselsInPort.add(new ArrayList<Integer>());
			for (int j = 0; j < vesselNum; j++) {
				if (ArrayUtils.contains(portOrder[j], i)) {
					vesselsInPort.get(i).add(j);
				}
			}
		}
		travelTimeVolocityLimit = new float[vesselNum][];
		for (int i = 0; i < vesselNum; i++) {
			travelTimeVolocityLimit[i] = new float[portOrder[i].length];
			for (int j = 0; j < portOrder[i].length; j++) {
				travelTimeVolocityLimit[i][j] = (travelTimeLimit[i][j][1] - travelTimeLimit[i][j][0]) * 0.15f;
			}
		}
		pc = new float[swarmNum];
		for (int i = 1; i <= swarmNum; i++) {
			pc[i - 1] = (float) (0.05 + 0.45 * ((Math.exp(10.0 * (i - 1) / (swarmNum - 1)) - 1) / (Math.exp(10) - 1)));
		}
	}

	public ArrayList<ArrayList<MyLinkedList<Integer>>> initBerthAllocation() {
		ArrayList<ArrayList<MyLinkedList<Integer>>> berthAllocation = new ArrayList<ArrayList<MyLinkedList<Integer>>>();
		for (int i = 0; i < portNum; i++) {
			berthAllocation.add(new ArrayList<MyLinkedList<Integer>>());
			ArrayList<ArrayList<Integer>> tempList = new ArrayList<ArrayList<Integer>>();
			for (int b = 0; b < berthNum[i]; b++) {
				berthAllocation.get(i).add(new MyLinkedList<Integer>());
				tempList.add(new ArrayList<Integer>());
			}
			TreeMap<Float, int[]> vesList = new TreeMap<>();
			for (int j = 0; j < vesselsInPort.get(i).size(); j++) {
				float min = 10000;
				float secMin = 10000;
				int idx = -1;
				for (int k = 0; k < berthNum[i]; k++) {
					if (min > handlingTimes[i][k][vesselsInPort.get(i).get(j)]) {
						secMin = min;
						min = handlingTimes[i][k][vesselsInPort.get(i).get(j)];
						idx = k;
					}
				}
				if (secMin - min >= 0.1 * secMin) {
					tempList.get(idx).add(vesselsInPort.get(i).get(j));
				} else {
					int ves = vesselsInPort.get(i).get(j);
					int p = ArrayUtils.indexOf(portOrder[ves], i);
					vesList.put(timeWindows[ves][p][1] - timeWindows[ves][p][0], new int[] { j, idx });
				}
			}
			while (vesList.size() > 0) {
				int v = vesList.get(vesList.lastKey())[0];
				int[] temp = new int[berthNum[i]];
				for (int j = 0; j < berthNum[i]; j++) {
					temp[j] = tempList.get(j).size();
				}
				int tempMin = Arrays.stream(temp).min().getAsInt();
				ArrayList<Integer> tempArr = new ArrayList<>();
				for (int j = 0; j < berthNum[i]; j++) {
					if (temp[j] == tempMin) {
						tempArr.add(j);
					}
				}
				tempArr.add(vesList.get(vesList.lastKey())[1]);
				double[] tempProb = new double[tempArr.size()];
				Arrays.fill(tempProb, 1.0 / tempArr.size());
				int b = new EnumeratedIntegerDistribution(tempArr.stream().mapToInt(Integer::valueOf).toArray(),
						tempProb).sample();
				tempList.get(b).add(vesselsInPort.get(i).get(v));
				vesList.remove(vesList.lastKey());
			}
			for (int j = 0; j < berthNum[i]; j++) {
				float[] temp = new float[tempList.get(j).size()];
				for (int k = 0; k < tempList.get(j).size(); k++) {
					int ves = tempList.get(j).get(k);
					int p = ArrayUtils.indexOf(portOrder[ves], i);
					temp[k] = (timeWindows[ves][p][1] + timeWindows[ves][p][0]) / 2;
				}
				float[] temp1 = ArrayUtils.clone(temp);
				Arrays.sort(temp1);
				for (int k = 0; k < tempList.get(j).size(); k++) {
					berthAllocation.get(i).get(j).add(tempList.get(j).get(ArrayUtils.indexOf(temp, temp1[k])));
				}
			}
		}
		berthAllocation = transferToFeasibleBerthAllocation(berthAllocation);
		return berthAllocation;
	}

	public ArrayList<ArrayList<MyLinkedList<Integer>>> transferToFeasibleBerthAllocation(
			ArrayList<ArrayList<MyLinkedList<Integer>>> berthAllocation) {
		ArrayList<ArrayList<MyLinkedList<Integer>>> berthAllocation1 = new ArrayList<ArrayList<MyLinkedList<Integer>>>();
		for (int i = 0; i < portNum; i++) {
			berthAllocation1.add(new ArrayList<MyLinkedList<Integer>>());
			for (int j = 0; j < berthNum[i]; j++) {
				berthAllocation1.get(i).add(new MyLinkedList<Integer>());
			}
		}
		int[] order = new int[vesselNum];
		int sum = 0;
		int sum1 = 0;
		int sum2 = 0;
		for (int i = 0; i < vesselNum; i++) {
			sum += portOrder[i].length;
		}
		while (sum1 < sum) {
			sum2 = sum1;
			sum1 = 0;
			for (int ves = 0; ves < vesselNum; ves++) {
				while (order[ves] < portOrder[ves].length) {
					int p = portOrder[ves][order[ves]];
					int b;
					for (b = 0; b < berthNum[p]; b++) {
						if (berthAllocation.get(p).get(b).size() != 0 && berthAllocation.get(p).get(b).get(0) == ves) {
							berthAllocation1.get(p).get(b).add(ves);
							berthAllocation.get(p).get(b).remove(0);
							order[ves]++;
							break;
						}
					}
					if (order[ves] == portOrder[ves].length || b == berthNum[p]) {
						break;
					}
				}
			}
			for (int i = 0; i < vesselNum; i++) {
				sum1 += order[i];
			}
			if (sum2 == sum1) {
				ArrayList<Integer> tempArr = new ArrayList<>();
				for (int i = 0; i < vesselNum; i++) {
					if (order[i] < portOrder[i].length) {
						tempArr.add(i);
					}
				}
				double[] tempProb = new double[tempArr.size()];
				Arrays.fill(tempProb, 1.0 / tempArr.size());
				int v = new EnumeratedIntegerDistribution(tempArr.stream().mapToInt(Integer::valueOf).toArray(),
						tempProb).sample();
				int p = portOrder[v][order[v]];
				for (int i = 0; i < berthNum[p]; i++) {
					if (berthAllocation.get(p).get(i).remove(Integer.valueOf(v))) {
						berthAllocation.get(p).get(i).addFirst(v);
						break;
					}
				}
			}
		}
		return berthAllocation1;
	}

	public float[][] initTravelTime() {
		UniformRealDistribution urd = new UniformRealDistribution(14, 16);
		float[][] travelTime = new float[vesselNum][];
		for (int i = 0; i < vesselNum; i++) {
			travelTime[i] = new float[portOrder[i].length];
			for (int j = 0; j < portOrder[i].length - 1; j++) {
				if (j == 0) {
					travelTime[i][j] = distances2[i] / (float) urd.sample();
				}
				travelTime[i][j + 1] = distances1[portOrder[i][j]][portOrder[i][j + 1]] / (float) urd.sample();
			}
		}
		return travelTime;
	}

	public float calCost(float[][] travelTime, ArrayList<ArrayList<MyLinkedList<Integer>>> berthAllocation) {
		float[][] arriveTime = new float[vesselNum][];
		float[][] handlingBeginTime = new float[vesselNum][];
		float[][] handlingEndTime = new float[vesselNum][];
		float[][] deltaEFTTime = new float[vesselNum][];
		int sum = 0;
		int sum1 = 0;
		int sum2 = 0;
		for (int i = 0; i < vesselNum; i++) {
			arriveTime[i] = new float[portOrder[i].length];
			arriveTime[i][0] = travelTime[i][0];
			handlingBeginTime[i] = new float[portOrder[i].length];
			handlingEndTime[i] = new float[portOrder[i].length];
			deltaEFTTime[i] = new float[portOrder[i].length];
			sum += portOrder[i].length;
		}
		float[][] berthTime = new float[portNum][];
		int[][] berthOrder = new int[portNum][];
		for (int i = 0; i < portNum; i++) {
			berthTime[i] = new float[berthNum[i]];
			berthOrder[i] = new int[berthNum[i]];
		}
		int[] order = new int[vesselNum];
		while (sum1 < sum) {
			sum2 = sum1;
			sum1 = 0;
			for (int ves = 0; ves < vesselNum; ves++) {
				while (order[ves] < portOrder[ves].length) {
					int p = portOrder[ves][order[ves]];
					int b;
					for (b = 0; b < berthNum[p]; b++) {
						if (berthAllocation.get(p).get(b).size() != berthOrder[p][b] && berthAllocation
								.get(portOrder[ves][order[ves]]).get(b).get(berthOrder[p][b]) == ves) {
							handlingBeginTime[ves][order[ves]] = (float) Arrays.stream(new double[] {
									timeWindows[ves][order[ves]][0], arriveTime[ves][order[ves]], berthTime[p][b] })
									.max().getAsDouble();
							handlingEndTime[ves][order[ves]] = handlingBeginTime[ves][order[ves]]
									+ handlingTimes[p][b][ves];
							deltaEFTTime[ves][order[ves]] = handlingEndTime[ves][order[ves]]
									- timeWindows[ves][order[ves]][1] < 0 ? 0
											: handlingEndTime[ves][order[ves]] - timeWindows[ves][order[ves]][1];
							if (order[ves] != portOrder[ves].length - 1) {
								arriveTime[ves][order[ves] + 1] = handlingEndTime[ves][order[ves]]
										+ travelTime[ves][order[ves] + 1];
							}
							berthTime[p][b] = handlingEndTime[ves][order[ves]];
							order[ves]++;
							berthOrder[p][b]++;
							break;
						}
					}
					if (order[ves] == portOrder[ves].length || b == berthNum[p]) {
						break;
					}
				}
			}
			for (int i = 0; i < vesselNum; i++) {
				sum1 += order[i];
			}
			if (sum1 == sum2) {
				return Float.MAX_VALUE;
			}
		}

		float idleCost = 0;
		float handlingCost = 0;
		float delayCost = 0;
		float fuelCost = 0;
		for (int i = 0; i < vesselNum; i++) {
			for (int j = 0; j < portOrder[i].length; j++) {
				idleCost += (handlingBeginTime[i][j] - arriveTime[i][j]);
				handlingCost += (handlingEndTime[i][j] - handlingBeginTime[i][j]);
				delayCost += deltaEFTTime[i][j];
				if (j == 0) {
					fuelCost += Fc * fd * Math.pow(distances2[i] / speedLimit[0], 3) * Math.pow(travelTime[i][j], -2);
				} else {
					fuelCost += Fc * fd * Math.pow(distances1[portOrder[i][j - 1]][portOrder[i][j]] / speedLimit[0], 3)
							* Math.pow(travelTime[i][j], -2);
				}
			}
		}
		idleCost *= Ic;
		handlingCost *= Hc;
		delayCost *= Dc;
//		System.out.println("idleCost:"+idleCost);
//		System.out.println("handlingCost:"+handlingCost);
//		System.out.println("delayCost"+delayCost);
//		System.out.println("fuelCost"+fuelCost);
		return idleCost + handlingCost + delayCost + fuelCost;
	}

	public int[][] twoOptType(int size) {
		int[][] tot = new int[size * (size - 1) / 2][2];
		int idx = 0;
		for (int i = 0; i < size - 1; i++) {
			for (int j = i + 2; j < size + 1; j++) {
				tot[idx][0] = i;
				tot[idx][1] = j;
				idx++;
			}
		}
		return tot;
	}

	public int[][] orOptType(int[] lens, int size) {
		int ootLen=0;
		for(int len:lens) {
			ootLen+=(size - len + 1) * (size - len);
		}
		int[][] oot = new int[ootLen][3];
		int idx = 0;
		for(int len:lens) {
			for (int i = 0; i < size - len + 1; i++) {
				for (int j = -1; j < i - 1; j++) {
					oot[idx][0]=len;
					oot[idx][1] = i;
					oot[idx][2] = j;
					idx++;
				}
				for (int j = i + len; j < size; j++) {
					oot[idx][0]=len;
					oot[idx][1] = i;
					oot[idx][2] = j;
					idx++;
				}
			}
		}
		return oot;
	}

	public void crossExchange(MyLinkedList<Integer> l1, MyLinkedList<Integer> l2, int[][] idxs) {
		Node<Integer> l1Beg = null;
		Node<Integer> l1End = null;
		Node<Integer> l1BegPPrev = null;
		Node<Integer> l1EndPNext = null;

		Node<Integer> l2Beg = null;
		Node<Integer> l2End = null;
		Node<Integer> l2BegPPrev = null;
		Node<Integer> l2EndPNext = null;

		if (idxs[0][1] == -1) {
			if (idxs[0][0] < l1.size) {
				l1EndPNext = l1.node(idxs[0][0]);
				l1BegPPrev = l1EndPNext.prev;
			} else {
				l1EndPNext = null;
				if (l1.size != 0) {
					l1BegPPrev = l1.node(l1.size - 1);
				}
			}
		} else {
			l1Beg = l1.node(idxs[0][0]);
			l1End = l1.node(idxs[0][1]);
			l1BegPPrev = l1Beg.prev;
			l1EndPNext = l1End.next;
		}

		if (idxs[1][1] == -1) {
			if (idxs[1][0] < l2.size) {
				l2EndPNext = l2.node(idxs[1][0]);
				l2BegPPrev = l2EndPNext.prev;
			} else {
				l2EndPNext = null;
				if (l2.size != 0) {
					l2BegPPrev = l2.node(l2.size - 1);
				}
			}
		} else {
			l2Beg = l2.node(idxs[1][0]);
			l2End = l2.node(idxs[1][1]);
			l2BegPPrev = l2Beg.prev;
			l2EndPNext = l2End.next;
		}

		if (l1BegPPrev == null && l1EndPNext == null) {
			l1.first = null;
			l1.last = null;
		}
		if (l2BegPPrev == null && l2EndPNext == null) {
			l2.first = null;
			l2.last = null;
		}

		if (l1Beg != null) {
			if (l2BegPPrev != null) {
				l2BegPPrev.next = l1Beg;
			} else {
				l2.first = l1Beg;
			}
			l1Beg.prev = l2BegPPrev;
		} else {
			if (l2BegPPrev != null) {
				l2BegPPrev.next = l2EndPNext;
			} else {
				l2.first = l2EndPNext;
			}
		}

		if (l1End != null) {
			if (l2EndPNext != null) {
				l2EndPNext.prev = l1End;
			} else {
				l2.last = l1End;
			}
			l1End.next = l2EndPNext;
		} else {
			if (l2EndPNext != null) {
				l2EndPNext.prev = l2BegPPrev;
			} else {
				l2.last = l2BegPPrev;
			}
		}

		if (l2Beg != null) {
			if (l1BegPPrev != null) {
				l1BegPPrev.next = l2Beg;
			} else {
				l1.first = l2Beg;
			}
			l2Beg.prev = l1BegPPrev;
		} else {
			if (l1BegPPrev != null) {
				l1BegPPrev.next = l1EndPNext;
			} else {
				l1.first = l1EndPNext;
			}
		}

		if (l2End != null) {
			if (l1EndPNext != null) {
				l1EndPNext.prev = l2End;
			} else {
				l1.last = l2End;
			}
			l2End.next = l1EndPNext;
		} else {
			if (l1EndPNext != null) {
				l1EndPNext.prev = l1BegPPrev;
			} else {
				l1.last = l1BegPPrev;
			}
		}
		for (int i = 0; i < 2; i++) {
			if (idxs[i][1] < 0) {
				idxs[i][0] = 0;
			}
		}
		int sizeChange = idxs[0][1] - idxs[0][0] - (idxs[1][1] - idxs[1][0]);
		l1.size -= sizeChange;
		l2.size += sizeChange;
	}

	public void icrossExchange(MyLinkedList<Integer> l1, MyLinkedList<Integer> l2, int[][] idxs) {
		Node<Integer> l1Beg = null;
		Node<Integer> l1End = null;
		Node<Integer> l1BegPPrev = null;
		Node<Integer> l1EndPNext = null;

		Node<Integer> l2Beg = null;
		Node<Integer> l2End = null;
		Node<Integer> l2BegPPrev = null;
		Node<Integer> l2EndPNext = null;

		if (idxs[0][1] == -1) {
			if (idxs[0][0] < l1.size) {
				l1EndPNext = l1.node(idxs[0][0]);
				l1BegPPrev = l1EndPNext.prev;
			} else {
				l1EndPNext = null;
				if (l1.size != 0) {
					l1BegPPrev = l1.node(l1.size - 1);
				}
			}
		} else {
			l1Beg = l1.node(idxs[0][0]);
			l1End = l1.node(idxs[0][1]);
			l1BegPPrev = l1Beg.prev;
			l1EndPNext = l1End.next;
			Node<Integer> temp = l1Beg;
			for (int i = 0; i < idxs[0][1] - idxs[0][0] + 1; i++) {
				Node<Integer> temp1 = temp.next;
				temp.next = temp.prev;
				temp.prev = temp1;
				temp = temp.prev;
			}
			temp = l1Beg;
			l1Beg = l1End;
			l1End = temp;
		}

		if (idxs[1][1] == -1) {
			if (idxs[1][0] < l2.size) {
				l2EndPNext = l2.node(idxs[1][0]);
				l2BegPPrev = l2EndPNext.prev;
			} else {
				l2EndPNext = null;
				if (l2.size != 0) {
					l2BegPPrev = l2.node(l2.size - 1);
				}
			}
		} else {
			l2Beg = l2.node(idxs[1][0]);
			l2End = l2.node(idxs[1][1]);
			l2BegPPrev = l2Beg.prev;
			l2EndPNext = l2End.next;
			Node<Integer> temp = l2Beg;
			for (int i = 0; i < idxs[1][1] - idxs[1][0] + 1; i++) {
				Node<Integer> temp1 = temp.next;
				temp.next = temp.prev;
				temp.prev = temp1;
				temp = temp.prev;
			}
			temp = l2Beg;
			l2Beg = l2End;
			l2End = temp;
		}

		if (l1BegPPrev == null && l1EndPNext == null) {
			l1.first = null;
			l1.last = null;
		}
		if (l2BegPPrev == null && l2EndPNext == null) {
			l2.first = null;
			l2.last = null;
		}

		if (l1Beg != null) {
			if (l2BegPPrev != null) {
				l2BegPPrev.next = l1Beg;
			} else {
				l2.first = l1Beg;
			}
			l1Beg.prev = l2BegPPrev;
		} else {
			if (l2BegPPrev != null) {
				l2BegPPrev.next = l2EndPNext;
			} else {
				l2.first = l2EndPNext;
			}
		}

		if (l1End != null) {
			if (l2EndPNext != null) {
				l2EndPNext.prev = l1End;
			} else {
				l2.last = l1End;
			}
			l1End.next = l2EndPNext;
		} else {
			if (l2EndPNext != null) {
				l2EndPNext.prev = l2BegPPrev;
			} else {
				l2.last = l2BegPPrev;
			}
		}

		if (l2Beg != null) {
			if (l1BegPPrev != null) {
				l1BegPPrev.next = l2Beg;
			} else {
				l1.first = l2Beg;
			}
			l2Beg.prev = l1BegPPrev;
		} else {
			if (l1BegPPrev != null) {
				l1BegPPrev.next = l1EndPNext;
			} else {
				l1.first = l1EndPNext;
			}
		}

		if (l2End != null) {
			if (l1EndPNext != null) {
				l1EndPNext.prev = l2End;
			} else {
				l1.last = l2End;
			}
			l2End.next = l1EndPNext;
		} else {
			if (l1EndPNext != null) {
				l1EndPNext.prev = l1BegPPrev;
			} else {
				l1.last = l1BegPPrev;
			}
		}
		for (int i = 0; i < 2; i++) {
			if (idxs[i][1] < 0) {
				idxs[i][0] = 0;
			}
		}
		int sizeChange = idxs[0][1] - idxs[0][0] - (idxs[1][1] - idxs[1][0]);
		l1.size -= sizeChange;
		l2.size += sizeChange;
	}

	public float[][] initTravelTimeVelocity() {
		float[][] velocity = new float[vesselNum][];
		for (int i = 0; i < vesselNum; i++) {
			velocity[i] = new float[portOrder[i].length];
			for (int j = 0; j < portOrder[i].length; j++) {
				velocity[i][j] = (float) new UniformRealDistribution(-travelTimeVolocityLimit[i][j],
						travelTimeVolocityLimit[i][j]).sample();
			}
		}
		return velocity;
	}

	public void getPBest() {
		if(PBest==null) {
			PBest=new float[swarmNum][][];
			for (int i = 0; i < swarmNum; i++) {
				PBestCost[i] = travelTimeSwarmCost[i];
				PBest[i] = new float[vesselNum][];
				for (int j = 0; j < vesselNum; j++) {
					PBest[i][j] = travelTimeSwarm[i][j].clone();
				}
			}
		}else {
			for (int i = 0; i < swarmNum; i++) {
				if (PBestCost[i] > travelTimeSwarmCost[i]) {
					flags[i]=0;
					PBestCost[i] = travelTimeSwarmCost[i];
					PBest[i] = new float[vesselNum][];
					for (int j = 0; j < vesselNum; j++) {
						PBest[i][j] = travelTimeSwarm[i][j].clone();
					}
				}else {
					flags[i]++;
				}
			}
		}
		
	}

	public void getGBest() {
		int idx=-1;
		for (int i = 0; i < swarmNum; i++) {
			if (travelTimeSwarmCost[i] < GBestCost) {
				GBestCost = travelTimeSwarmCost[i];
				idx=i;
			}
		}
		if(idx!=-1) {
			GBest = new float[vesselNum][];
			for (int j = 0; j < vesselNum; j++) {
				GBest[j] = travelTimeSwarm[idx][j].clone();
			}
		}
	}
	
	public void getBest() {
		if(GBestCost<bestCost || VNSBestCost<bestCost) {
			bestCost=GBestCost;
			bestTravelTime=new float[vesselNum][];
			for (int j = 0; j < vesselNum; j++) {
				bestTravelTime[j] = GBest[j].clone();
			}
			bestBerthAllocation=new ArrayList<ArrayList<MyLinkedList<Integer>>>();
			for(int i=0;i<portNum;i++) {
				bestBerthAllocation.add(new ArrayList<MyLinkedList<Integer>>());
				for(int j=0;j<berthNum[i];j++) {
					bestBerthAllocation.get(i).add(VNSBestBerthAllocation.get(i).get(j).clone());
				}
			}
		}
	}

	public void CLPSOPositionUpdating(float omega, int idx, Random rand) {
		float pc = rand.nextFloat();
		float[][] tempTravelTime = new float[vesselNum][];
//		int[] choose=new int[swarmNum];
//		for(int i=0;i<swarmNum;i++) {
//			choose[i]=i;
//		}
//		ArrayUtils.shuffle(choose);
		int idxPBest = -1;
		if (pc > this.pc[idx]) {
			idxPBest = idx;
		} else {
			int vel1 = rand.nextInt(swarmNum);
			int vel2 = rand.nextInt(swarmNum);
			while (vel1 == vel2 || vel1 == idx || vel2 == idx) {
				vel1 = rand.nextInt(swarmNum);
				vel2 = rand.nextInt(swarmNum);
			}
			idxPBest = PBestCost[vel1] > PBestCost[vel2] ? vel2 : vel1;
		}
		boolean flag = true;
		t: for (int v = 0; v < vesselNum; v++) {
			tempTravelTime[v] = new float[portOrder[v].length];
			for (int p = 0; p < portOrder[v].length; p++) {
				travelTimeSwarmVelocity[idx][v][p] = omega * travelTimeSwarmVelocity[idx][v][p]
						+ c1 * rand.nextFloat() * (PBest[idx][v][p] - travelTimeSwarm[idxPBest][v][p]);
				if (travelTimeSwarmVelocity[idx][v][p] > travelTimeVolocityLimit[v][p]) {
					travelTimeSwarmVelocity[idx][v][p] = travelTimeVolocityLimit[v][p];
				} else if (travelTimeSwarmVelocity[idx][v][p] < -travelTimeVolocityLimit[v][p]) {
					travelTimeSwarmVelocity[idx][v][p] = -travelTimeVolocityLimit[v][p];
				}
				float temp = travelTimeSwarm[idx][v][p] + travelTimeSwarmVelocity[idx][v][p];
				if (temp > travelTimeLimit[v][p][0] && temp < travelTimeLimit[v][p][1]) {
					tempTravelTime[v][p] = temp;
				} else {
					flag = false;
					break t;
				}
			}
		}
		if (flag) {
			travelTimeSwarm[idx] = tempTravelTime;
			travelTimeSwarmCost[idx] = calCost(tempTravelTime, VNSBestBerthAllocation);
		}
	}

	public void PSOPositionUpdating(float omega, int idx, Random rand) {
		float[][] tempTravelTime = new float[vesselNum][];
		boolean flag = true;
		t: for (int v = 0; v < vesselNum; v++) {
			tempTravelTime[v] = new float[portOrder[v].length];
			for (int p = 0; p < portOrder[v].length; p++) {
				travelTimeSwarmVelocity[idx][v][p] = omega * travelTimeSwarmVelocity[idx][v][p]
						+ c1 * rand.nextFloat() * (PBest[idx][v][p] - travelTimeSwarm[idx][v][p])
						+ c2 * rand.nextFloat() * (GBest[v][p] - travelTimeSwarm[idx][v][p]);
				if (travelTimeSwarmVelocity[idx][v][p] > travelTimeVolocityLimit[v][p]) {
					travelTimeSwarmVelocity[idx][v][p] = travelTimeVolocityLimit[v][p];
				} else if (travelTimeSwarmVelocity[idx][v][p] < -travelTimeVolocityLimit[v][p]) {
					travelTimeSwarmVelocity[idx][v][p] = -travelTimeVolocityLimit[v][p];
				}
				float temp = travelTimeSwarm[idx][v][p] + travelTimeSwarmVelocity[idx][v][p];
				if (temp > travelTimeLimit[v][p][0] && temp < travelTimeLimit[v][p][1]) {
					tempTravelTime[v][p] = temp;
				} else {
					flag = false;
					break t;
				}
			}
		}
		if (flag) {
			travelTimeSwarm[idx] = tempTravelTime;
			travelTimeSwarmCost[idx] = calCost(tempTravelTime, VNSBestBerthAllocation);
		}
	}

	public ImmutablePair<ArrayList<ArrayList<MyLinkedList<Integer>>>, int[][]> shaking(int k) {
		ArrayList<ArrayList<MyLinkedList<Integer>>> berthAllocation = new ArrayList<ArrayList<MyLinkedList<Integer>>>();
		k = k % 6;
		int[][] choosedBerth = new int[portNum][2];
		for (int i = 0; i < portNum; i++) {
			int[] tempBerth = new int[berthNum[i]];
			berthAllocation.add(new ArrayList<MyLinkedList<Integer>>());
			for (int j = 0; j < berthNum[i]; j++) {
				tempBerth[j] = j;
				berthAllocation.get(i).add(this.berthAllocation.get(i).get(j).clone());
			}
			ArrayUtils.shuffle(tempBerth);
			choosedBerth[i][0] = tempBerth[0];
			choosedBerth[i][1] = tempBerth[1];
			int subLen1;
			int subLen2;
			switch (k) {
			case 0:
				subLen1 = Arrays.stream(new int[] { berthAllocation.get(i).get(tempBerth[0]).size(), 1 }).min()
						.getAsInt();
				subLen1 = new UniformIntegerDistribution(0, subLen1).sample();
				subLen2 = Arrays.stream(new int[] { berthAllocation.get(i).get(tempBerth[1]).size(), 1 }).min()
						.getAsInt();
				subLen2 = new UniformIntegerDistribution(0, subLen2).sample();
				break;
			case 1:
				subLen1 = Arrays.stream(new int[] { berthAllocation.get(i).get(tempBerth[0]).size(), 2 }).min()
						.getAsInt();
				subLen1 = new UniformIntegerDistribution(0, subLen1).sample();
				subLen2 = Arrays.stream(new int[] { berthAllocation.get(i).get(tempBerth[1]).size(), 2 }).min()
						.getAsInt();
				subLen2 = new UniformIntegerDistribution(0, subLen2).sample();
				break;
			case 2:
				subLen1 = Arrays.stream(new int[] { berthAllocation.get(i).get(tempBerth[0]).size(), 3 }).min()
						.getAsInt();
				subLen1 = new UniformIntegerDistribution(0, subLen1).sample();
				subLen2 = Arrays.stream(new int[] { berthAllocation.get(i).get(tempBerth[1]).size(), 3 }).min()
						.getAsInt();
				subLen2 = new UniformIntegerDistribution(0, subLen2).sample();
				break;
			case 3:
				subLen1 = Arrays.stream(new int[] { berthAllocation.get(i).get(tempBerth[0]).size(), 4 }).min()
						.getAsInt();
				subLen1 = new UniformIntegerDistribution(0, subLen1).sample();
				subLen2 = Arrays.stream(new int[] { berthAllocation.get(i).get(tempBerth[1]).size(), 4 }).min()
						.getAsInt();
				subLen2 = new UniformIntegerDistribution(0, subLen2).sample();
				break;
			case 4:
				subLen1 = Arrays.stream(new int[] { berthAllocation.get(i).get(tempBerth[0]).size(), 5 }).min()
						.getAsInt();
				subLen1 = new UniformIntegerDistribution(0, subLen1).sample();
				subLen2 = Arrays.stream(new int[] { berthAllocation.get(i).get(tempBerth[1]).size(), 5 }).min()
						.getAsInt();
				subLen2 = new UniformIntegerDistribution(0, subLen2).sample();
				break;
			case 5:
				subLen1 = berthAllocation.get(i).get(tempBerth[0]).size();
				subLen1 = new UniformIntegerDistribution(0, subLen1).sample();
				subLen2 = berthAllocation.get(i).get(tempBerth[1]).size();
				subLen2 = new UniformIntegerDistribution(0, subLen2).sample();
				break;
			default:
				subLen1 = -1;
				subLen2 = -1;
				break;
			}
			int beg1 = new UniformIntegerDistribution(0, berthAllocation.get(i).get(tempBerth[0]).size() - subLen1)
					.sample();
			int beg2 = new UniformIntegerDistribution(0, berthAllocation.get(i).get(tempBerth[1]).size() - subLen2)
					.sample();
			float icrossPosiblility = 1.0f / vesselsInPort.get(i).size();
			float rand = new Random().nextFloat();
			if (subLen1 == 0) {
				subLen1 = -beg1;
			}
			if (subLen2 == 0) {
				subLen2 = -beg2;
			}
			if (rand > icrossPosiblility) {
				crossExchange(berthAllocation.get(i).get(tempBerth[0]), berthAllocation.get(i).get(tempBerth[1]),
						new int[][] { new int[] { beg1, beg1 + subLen1 - 1 }, new int[] { beg2, beg2 + subLen2 - 1 } });
			} else {
				icrossExchange(berthAllocation.get(i).get(tempBerth[0]), berthAllocation.get(i).get(tempBerth[1]),
						new int[][] { new int[] { beg1, beg1 + subLen1 - 1 }, new int[] { beg2, beg2 + subLen2 - 1 } });
			}
		}
		ImmutablePair<ArrayList<ArrayList<MyLinkedList<Integer>>>, int[][]> ip = new ImmutablePair<ArrayList<ArrayList<MyLinkedList<Integer>>>, int[][]>(
				berthAllocation, choosedBerth);
		return ip;
	}

	public ImmutablePair<ArrayList<ArrayList<MyLinkedList<Integer>>>, Float> localSearch(ImmutablePair<ArrayList<ArrayList<MyLinkedList<Integer>>>, int[][]> ip) {
		ArrayList<ArrayList<MyLinkedList<Integer>>> berthAllocation=ip.left;
		int[][] choosedBerth=ip.right;
		ArrayList<ArrayList<MyLinkedList<Integer>>> choosedCopy=new ArrayList<ArrayList<MyLinkedList<Integer>>>();
		int[] portOrder=new int[portNum];
		for(int i=0;i<portNum;i++) {
			portOrder[i]=i;
			choosedCopy.add(new ArrayList<MyLinkedList<Integer>>());
			for(int j=0;j<2;j++) {
				choosedCopy.get(i).add(berthAllocation.get(i).get(choosedBerth[i][j]).clone());
			}
		}
		ArrayUtils.shuffle(portOrder);
		int[][] orOptType=null;
		int[][] twoOptType=null;
		double localSearchCost=Float.MAX_VALUE;
		double localSearchBestCost=Float.MAX_VALUE;
		for(int i:portOrder) {
			for(int j=0;j<2;j++) {
				MyLinkedList<Integer> localSearch=choosedCopy.get(i).get(j);
				int localSearchLen=localSearch.size();
				MyLinkedList<Integer> localSearchBest=null;
				localSearchCost=calCost(GBest, berthAllocation);
				if(localSearchCost<=localSearchBestCost) {
					localSearchBestCost=localSearchCost;
					localSearchBest=localSearch;
				}
				float rand=new Random().nextFloat();
				if(rand<0.5) {
					if(localSearchLen==2) {
						orOptType=orOptType(new int[] {1}, localSearchLen);
					}else if(localSearchLen==3) {
						orOptType=orOptType(new int[] {1,2}, localSearchLen);
					}else if(localSearchLen>3) {
						orOptType=orOptType(new int[] {1,2,3}, localSearchLen);
					}
					if(orOptType!=null) {
						for(int k=0;k<orOptType.length;k++) {
							localSearch=choosedCopy.get(i).get(j).clone();
							localSearch.orOpt(orOptType[k][0], orOptType[k][1], orOptType[k][2]);
							berthAllocation.get(i).set(choosedBerth[i][j], localSearch);
							localSearchCost=calCost(GBest, berthAllocation);
							if(localSearchCost<localSearchBestCost) {
								localSearchBestCost=localSearchCost;
								localSearchBest=localSearch;
							}
						}
						orOptType=null;
					}
				}else {
					if(localSearchLen>1) {
						twoOptType=twoOptType(localSearchLen);
					}
					if(twoOptType!=null) {
						for(int k=0;k<twoOptType.length;k++) {
							localSearch=choosedCopy.get(i).get(j).clone();
							localSearch.twoOpt(twoOptType[k][0],twoOptType[k][1]);
							berthAllocation.get(i).set(choosedBerth[i][j], localSearch);
							localSearchCost=calCost(GBest, berthAllocation);
							if(localSearchCost<localSearchBestCost) {
								localSearchBestCost=localSearchCost;
								localSearchBest=localSearch;
							}
						}
						twoOptType=null;
					}
				}
				berthAllocation.get(i).set(choosedBerth[i][j], localSearchBest);
			}
		}
		ImmutablePair<ArrayList<ArrayList<MyLinkedList<Integer>>>, Float> ip1=new ImmutablePair<>(berthAllocation, (float)localSearchBestCost);
		return ip1;
	}
	
	public int updateBerthAllocation(int k,int iteration,int iterationNum) {
		ImmutablePair<ArrayList<ArrayList<MyLinkedList<Integer>>>, int[][]> ip=shaking(k);
		ImmutablePair<ArrayList<ArrayList<MyLinkedList<Integer>>>, Float> ip1=localSearch(ip);
		ArrayList<ArrayList<MyLinkedList<Integer>>> berthAllocation=ip1.left;
		float localSearchBestCost=ip1.right;
		float T0=5;
		float T=T0-T0*(iteration/1000)*1000/iterationNum;
		if(localSearchBestCost<VNSLastCost) {
			VNSLastCost=localSearchBestCost;
			this.berthAllocation=berthAllocation;
			if(localSearchBestCost<VNSBestCost) {
				VNSBestCost=localSearchBestCost;
				VNSBestBerthAllocation=new ArrayList<ArrayList<MyLinkedList<Integer>>>();
				for(int i=0;i<portNum;i++) {
					VNSBestBerthAllocation.add(new ArrayList<MyLinkedList<Integer>>());
					for(int j=0;j<berthNum[i];j++) {
						VNSBestBerthAllocation.get(i).add(berthAllocation.get(i).get(j).clone());
					}
				}
			}
			k=0;
		}else {
			float rand=new Random().nextFloat();
			if(rand<=Math.exp((VNSLastCost-localSearchBestCost)/T)) {
				VNSLastCost=localSearchBestCost;
				this.berthAllocation=berthAllocation;
				k=0;
			}else {
				k++;
			}
		}
		return k;
	}
	
	

	public void run() {
		travelTimeSwarm = new float[swarmNum][][];
		travelTimeSwarmCost = new float[swarmNum];
		travelTimeSwarmVelocity = new float[swarmNum][][];
		for (int i = 0; i < swarmNum; i++) {
			travelTimeSwarm[i] = initTravelTime();
			travelTimeSwarmVelocity[i] = initTravelTimeVelocity();

		}
		berthAllocation = initBerthAllocation();
		VNSBestBerthAllocation=new ArrayList<ArrayList<MyLinkedList<Integer>>>();
		for(int i=0;i<portNum;i++) {
			VNSBestBerthAllocation.add(new ArrayList<MyLinkedList<Integer>>());
			for(int j=0;j<berthNum[i];j++) {
				VNSBestBerthAllocation.get(i).add(berthAllocation.get(i).get(j).clone());
			}
		}
		PBestCost = new float[swarmNum];
		PBest = null;
		GBestCost = Float.MAX_VALUE;
		GBest = null;
		VNSLastCost=Float.MAX_VALUE;
		VNSBestCost=Float.MAX_VALUE;
		bestCost=Float.MAX_VALUE;
		flags = new int[swarmNum];
		PSOiterationNum=2000;
		VNSiterationNum=1000;
		Random rand = new Random();
		for (int i = 0; i < iterationNum; i++) {
			for(int s=0;s<swarmNum;s++) {
				travelTimeSwarmCost[s]=calCost(travelTimeSwarm[s], VNSBestBerthAllocation);
			}
			for(int j=0;j<PSOiterationNum;j++) {
				getPBest();
				getGBest();
				float omega = (float) (0.9 - 0.5 * i / PSOiterationNum);
				for (int k = 0; k < swarmNum; k++) {
					if (flags[k] < rg) {
						CLPSOPositionUpdating(omega, k, rand);
					} else {
						PSOPositionUpdating(omega, k, rand);
						flags[k] = 0;
					}
				}
//				System.out.println("PSOiterationNum "+j+" "+GBestCost);
			}
			System.out.println("iterationNum "+i+"PSOiterationNum "+" "+GBestCost);
			System.out.println(Arrays.deepToString(GBest));
			System.out.println(VNSBestBerthAllocation.toString());
			getBest();
			int k=0;
			for(int j=0;j<VNSiterationNum;j++) {
				k=updateBerthAllocation(k, j, VNSiterationNum);
			}
			System.out.println("iterationNum "+i+"VNSiterationNum "+" "+VNSBestCost);
			getBest();
//			VNSLastCost=M;
//			VNSBestCost=M;
//			PBest=null;
//			GBest=null;
//			GBestCost=M;
		}
		System.out.println(bestCost);
	}
}
