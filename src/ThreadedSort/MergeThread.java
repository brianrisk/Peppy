package ThreadedSort;

import java.util.ArrayList;
import java.util.List;


public class MergeThread<E extends Comparable<E>> extends Thread {
	
	List<E> result;
	List<E> listA;
	List<E> listB;
	
	public MergeThread(List<E> listA, List<E> listB) {
		this.listA = listA;
		this.listB = listB;
	}

	public void run() {
		int sizeA = listA.size();
		int sizeB = listB.size();
		int indexA = 0;
		int indexB = 0;
		result = new ArrayList<E>(sizeA + sizeB);
		while (indexA < sizeA && indexB < sizeB) {
			if (listA.get(indexA).compareTo(listB.get(indexB)) <= 0) {
				result.add(listA.get(indexA));
				indexA++;
			} else {
				result.add(listB.get(indexB));
				indexB++;
			}
		}
		if (indexA < sizeA) {
			for (int i = indexA; i < sizeA; i++) {
				result.add(listA.get(i));
			}
		} else {
			for (int i = indexB; i < sizeB; i++) {
				result.add(listB.get(i));
			}
		}	
	}

	/**
	 * @return the ourList
	 */
	public List<E> getList() {
		return result;
	}

}