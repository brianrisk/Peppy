package ThreadedSort;

import java.util.*;

public class Sorter extends Thread {
	
	private static final int INSERTIONSORT_THRESHOLD = 7;
	
	private Object[] src;
	private Object[] dest;
	private int start;
	private int stop;
	
	/**
	 * Tester
	 */
//	public static void main(String args[]) {
//		int ourListSize = 27957796;
////		Object [] cat = null;
////		Arrays.sort(cat);
////		int ourListSize = 63;
//		Random random = new Random();
//		
//		//initialize the array
//		ArrayList<Integer> ourList = new ArrayList<Integer>(ourListSize);
//		for (int i = 0; i < ourListSize; i++) {
//			ourList.add(new Integer(random.nextInt(ourListSize)));
//		}
//
//		//time our sort
//		long startTime = System.currentTimeMillis();
//		sort(ourList, 16);
//		System.out.println("here's the time:  " + (System.currentTimeMillis() - startTime));
//		
////		for (Integer dude: ourList) {
////			System.out.println(dude);
////		}
//		
//		//Verify that we are sorting properly
//		int val= 0;
//		Integer dude;
//		for (int i = 0; i < ourListSize; i++) {
//			dude  = ourList.get(i);
//			if (dude.intValue() < val) {System.out.println("wrong: " + i + ", value: " + dude.intValue()); break;}
//			val = dude.intValue();
//		}
//		
//		//re-initialize the array
//		ourList = new ArrayList<Integer>(ourListSize);
//		for (int i = 0; i < ourListSize; i++) {
//			ourList.add(new Integer(random.nextInt(ourListSize)));
//		}
//		
//		//time a standard sort
//		startTime = System.currentTimeMillis();
//		Collections.sort(ourList);
//		System.out.println("here's the time:  " + (System.currentTimeMillis() - startTime));
//	}
	
	
	public Sorter(Object[] src,
			Object[] dest,
			int start,
			int stop) {
		this.src = src;
		this.dest = dest;
		this.start = start;
		this.stop = stop;
//		System.out.println("new sort object.  start: " + start + ". stop: " + stop);
	}
	
	//@SuppressWarnings("unchecked")
	public void run() {
		mergeSort(src, dest, start, stop);
	}

	/**
	 * Based on java.util.Collections
	 * @param <T>
	 * @param list
	 * @param threadCount This should be a power of 2.  If it is not, 
	 * the parameter used will be the greatest power of 2 
	 * that is less than this number.
	 */
	public static <T extends Comparable<? super T>> void sort(List<T> list, int threadCount) {
		//make sure threadCount is power of 2
		double exponent = Math.log(threadCount) / Math.log(2);
		exponent = Math.floor(exponent);
		threadCount = (int) Math.pow(2, exponent);
		
		//convert list to array.  clone array.
		Object[] dest = list.toArray();
		Object[] src = (Object[])dest.clone();
		
		//This holds all of our threads which will be sorting
		ArrayList<Sorter> sortThreads = new ArrayList<Sorter>(threadCount);
		
		int start = 0;
		int stop = 0;
		int subUnitSize = list.size() / threadCount;
		for (int i = 0; i < threadCount - 1; i++) {
			start = i * subUnitSize;
			stop = (i + 1) * subUnitSize;

			/*
			 * Now that a sub unit of ourList is defined, we add it to a sorting thread,
			 * start that thread, and then add the running thread to a group of threads.
			 */
			Sorter thread = new Sorter(src, dest, start, stop);
			thread.start();
			sortThreads.add(thread);
		}
		
		// add the final thread which goes all the way to the end of the array
		Sorter thread = new Sorter(src, dest, stop, dest.length);
		thread.start();
		sortThreads.add(thread);
		
		//Wait for all threads to finish executing
		for (Sorter sorter: sortThreads) {
			try {
				sorter.join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		
		//Merge our sorted sections
		int mid;
		int groupCount = threadCount / 2;
		while (groupCount > 0) {
			/*
			 * Here we're grouping sections leaving the last, oddly shaped
			 * group to be picked up and merged after this for loop
			 */
			int lastGroup = groupCount;
			for (int groupNumber = 0; groupNumber < lastGroup; groupNumber++) {
				start = (groupNumber * 2 + 0) * subUnitSize;
				mid   = (groupNumber * 2 + 1) * subUnitSize;
				stop  = (groupNumber * 2 + 2) * subUnitSize;
				merge(src, dest, start, mid, stop);
			}
			start = lastGroup * subUnitSize;
			mid   = (lastGroup + 1) * subUnitSize;
			stop  = dest.length;
			merge(dest, src, start, mid, stop);
			
			groupCount /= 2;
			subUnitSize *= 2;
		}
		

		//set final order of list
		ListIterator<T> iterator = list.listIterator();
		for (int j=0; j<src.length; j++) {
			iterator.next();
			iterator.set((T)src[j]);
		}
	}
	
	private static void merge(
			Object[] dest,
			Object[] src,
			int start,
			int mid,
			int stop) {
		for(int i = start, p = start, q = mid; i < stop; i++) {
			if (q >= stop || p < mid && ((Comparable)src[p]).compareTo(src[q])<=0)
				dest[i] = src[p++];
			else
				dest[i] = src[q++];
		}
	}

	
//	public static void sort(Object[] a) {
//		Object[] aux = (Object[])a.clone();
//		mergeSort(aux, a, 0, a.length);
//	}

	
	/**
	 * Based on java.lang.Arrays
	 */
	private void mergeSort(Object[] src,
			Object[] dest,
			int start,
			int stop) {
		int length = stop - start;

		// Insertion sort on smallest arrays
		if (length < INSERTIONSORT_THRESHOLD) {
			for (int i=start; i<stop; i++)
				for (int j=i; j>start &&
				((Comparable) dest[j-1]).compareTo(dest[j])>0; j--)
					swap(dest, j, j-1);
			return;
		}

		// Recursively sort halves of dest into src
		int mid = (start + stop) >>> 1;
		mergeSort(dest, src, start, mid);
		mergeSort(dest, src, mid, stop);

		// If list is already sorted, just copy from src to dest.  This is an
		// optimization that results in faster sorts for nearly ordered lists.
		if (((Comparable)src[mid-1]).compareTo(src[mid]) <= 0) {
			System.arraycopy(src, start, dest, start, length);
			return;
		}

		// Merge sorted halves (now in src) into dest
		for(int i = start, p = start, q = mid; i < stop; i++) {
			if (q >= stop || p < mid && ((Comparable)src[p]).compareTo(src[q])<=0)
				dest[i] = src[p++];
			else
				dest[i] = src[q++];
		}
	}

	
	private static void swap(Object[] x, int a, int b) {
		Object t = x[a];
		x[a] = x[b];
		x[b] = t;
	}

}
