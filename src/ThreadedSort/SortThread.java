package ThreadedSort;

import java.util.Collections;
import java.util.List;


/*
 * a Runnable object which is initialized with a collection,
 * sorts that collection in its thread.  returns that
 * sorted collection when queried.
 */
public class SortThread<E extends Comparable<E>> extends Thread {
	
	List<E> list;
	
	public SortThread(List<E> list) {
		this.list = list;
	}

	//@SuppressWarnings("unchecked")
	public void run() {
		Collections.sort(list);
	}

	/**
	 * @return the ourList
	 */
	public List<E> getList() {
		return list;
	}

}
