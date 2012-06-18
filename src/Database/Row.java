package Database;

import java.io.File;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;

import Peppy.U;


/**
 * This is a generic data object.
 * 
 * There is a set of predefine columns.  However,
 * a row with any key value can potentially be added.
 * 
 * Copyright 2012, Brian Risk
 * Released under the Netscape Public License
 * 
 * @author Brian Risk
 *
 */
public class Row {
	
	protected static Hashtable<String, Class<?>> columns = new Hashtable<String, Class<?>>();
	protected Hashtable<String, Object> values = new Hashtable<String, Object>();


	public void set(String name, Object value) {
		/* okay, special case for peptides */
		if (name.equals("peptideSequence")) {
			String peptideSequence = (String) value;
			/* removing the stop so that it will match with peptides from protein sequences */
			if (peptideSequence.endsWith(".")) {
				peptideSequence = peptideSequence.substring(0, peptideSequence.length() - 1);
			}
			value = peptideSequence;
		}
		
		name = name.toLowerCase();
		Class<?> propertyType = columns.get(name);
		if (propertyType == null) {
			values.put(name, value);
		} else {
			if (propertyType.isInstance(value)) {
				values.put(name, value);
			} else {
				U.p("Value is not of type expected.  Expected:  " + propertyType.toString() + " Received: " + value.getClass());
				System.exit(1);
			}
		}
	}
	
	public Object get(String name) {
		Object out = values.get(name.toLowerCase());
		if (out == null) {
			
		}
		return out;
	}
	
	public File getFile(String name) {
		return (File) get(name);
	}
	
	public String getString(String name) {
		String out =  (String) get(name);
		return out;
	}
	
	public int getInt(String name) {
		return (Integer) get(name);
	}
	
	public double getDouble(String name) {
		return (Double) get(name);
	}
	
	public boolean getBoolean(String name) {
		return (Boolean) get(name);
	}
	
	public Hashtable<String, Class<?>> getColumns() {
		return columns;
	}
	
	
	public static Enumeration<String> getColumnNames() {
		return columns.keys();
	}

	

	
	public String toString() {
		StringBuffer out = new StringBuffer();
		Enumeration<Object> e = values.elements();
		Object nextElement;
		while (e.hasMoreElements()) {
			nextElement = e.nextElement();
			out.append(nextElement);
			if(e.hasMoreElements()) {
				out.append('\t');
			}
		}
		return out.toString();
	}
	
	
	

	
	/**
	 * Allows a user to define which fields are compared for equality
	 * @param other
	 * @param keys
	 * @return
	 */
	public boolean equals(Row other, ArrayList<String> keys) {
		for (String key: keys) {
			if (!other.getString(key).equals(getString(key))) return false;
		}
		return true;
	}
	
	
	/**
	 * Allows a user to define which field is compared to determine equality
	 * @param other
	 * @param keys
	 * @return
	 */
	public boolean equals(Row other, String key) {
		return other.getString(key).equals(getString(key));
	}
	
	public int getSize() {
		return values.size();
	}
	
	
}
