package Database;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;


/**
 * A take on the Table database concept wherein multiple rows can be packed
 * into one key value.
 * 
 * This allows "overriding" for the key value.  For example, the key value
 * could be "name"; if there happened to be two "Brian"s, both would be stored
 * for the key value of "Brian".
 * 
 * Copyright 2012, Brian Risk
 * 
 * 
 * @author Brian Risk
 *
 */
public class Table { 
	
	/* the name of our table */
	private String name = "";
	
	/* the predefined columns for this table */
	public static Hashtable<String, Column> columns = new Hashtable<String, Column>();
	
	/* where we store all the matches */
	Hashtable<String, ArrayList<Row>> rows = new Hashtable<String, ArrayList<Row>>();
	
	/*
	 * In some cases, when we add new rows, we want to keep only
	 * the "best" of them for any given key.  We need to set this
	 * property as well as define the column by which the best is
	 * measured.
	 */
	private boolean keepOnlyTheBest = false;
	private Column bestColumn;
	
	/* every table has a "key" value.  This is the column that must contain unique values.
	 * In fact, it is the key for our "rows" variable;  That is, if we were to cal
	 */
	private Column key;
	
	
	public Table(Column key) {
		this.key = key;
	}
	
	public Table(String name, Column key) {
		this.name = name;
		this.key = key;
	}
	

	
	
	/**
	 * Handles all of the headache of adding to a hashtable of array lists.
	 * @param key
	 * @param row
	 */
	public void add(Row row) {
		String keyValue = row.getString(key.getName());
		ArrayList<Row> existingMatches = get(keyValue);
		if (existingMatches == null) {
			ArrayList<Row> matchArray = new ArrayList<Row>();
			matchArray.add(row);
			rows.put(keyValue, matchArray);
		} else {
			
			/*
			 * if we are keeping only the best then we must see what the score
			 * of the existing matches are (and they should all be the same if
			 * there are more than one).  
			 */
			if (keepOnlyTheBest) {
				
				boolean isBetter = false;
				boolean isEqual = false;
				
				/* if our "best" feature is an integer */
				if (bestColumn.getType().equals(Integer.class)) {
					if (row.getInt(bestColumn.getName()) > existingMatches.get(0).getInt(bestColumn.getName())) isBetter = true;
					if (row.getInt(bestColumn.getName()) == existingMatches.get(0).getInt(bestColumn.getName())) isEqual = true;
				}
				
				/* if it is a double */
				if (bestColumn.getType().equals(Double.class)) {
					if (row.getDouble(bestColumn.getName()) > existingMatches.get(0).getDouble(bestColumn.getName())) isBetter = true;
					if (row.getDouble(bestColumn.getName()) == existingMatches.get(0).getDouble(bestColumn.getName())) isEqual = true;
				}
				
				/* If the new match has a higher score then we create a new ArrayList and usurp the existing with it. */
				if (isBetter) {
					ArrayList<Row> matchArray = new ArrayList<Row>();
					matchArray.add(row);
					rows.put(key.getName(), matchArray);
				}
				
				/* If the new score equals the existing score, we add it to the list */
				if (isEqual) {
					existingMatches.add(row);
				}
				
			} else {
				existingMatches.add(row);
			}
		}
	}
	
	
	/**
	 * 
	 */
	public ArrayList<Row> get(String key) {
		return rows.get(key);
	}
	
	public ArrayList<Row> get(Column key) {
		return get(key.getName());
	}
	
	
	/**
	 * 
	 * @return returns a full, unpacked ArrayList of all of our rows
	 */
	public ArrayList<Row> getRows() {
		/* determining the size of the output */
		int size = 0;
		Enumeration<ArrayList<Row>> e = rows.elements();
		while (e.hasMoreElements()) {
			size += e.nextElement().size();
		}
		
		/* create the output */
		ArrayList<Row> out = new ArrayList<Row>(size);
		e = rows.elements();
		while (e.hasMoreElements()) {
			out.addAll(e.nextElement());
		}
		
		return out;
	}
	
	public Column getKey() {
		return key;
	}
	
	
	public String getName() {
		return name;
	}




	/**
	 * 
	 * @return an ArrayList of all of the keys
	 */
	public ArrayList<String> getKeyValues() {
		return new ArrayList<String>(rows.keySet());
	}

	
	
	/**
	 * Writes our table to a file. using the predefined columns
	 * @param file
	 */
	public void write(File file) {
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
			Enumeration<Column> columnElements = columns.elements();
			
			ArrayList<Class<?>> types = new ArrayList<Class<?>>();
			ArrayList<String> names = new ArrayList<String>();
			Column column;
			while (columnElements.hasMoreElements()) {
				column = columnElements.nextElement();
				types.add(column.getType());
				names.add(column.getName());
			}
			
			/* our line we'll build and write */
			StringBuffer line;
			
			/* print the column types */
			line = new StringBuffer();
			for (Class<?> type: types) {
				line.append(type.getName());
				line.append('\t');
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	/**
	 * Loads a file of data into this table.  The file, at the very least
	 * must contain a column that matches in both name and data type to our
	 * key column
	 * 
	 * If a column matches that of a column of this table, then all
	 * of the data in that column must conform to the defined properties.
	 * 
	 * @param file
	 */
	public void load(File file) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			
			/* read the header lines */
			br.readLine();
			br.readLine();
			br.readLine();
			
			/* find out what property each column represents */
			String line = br.readLine();
			String [] propertyNames = line.split("\t");
			
			/* read in the first line */
			line = br.readLine();
			
			while (line != null) {
				String [] chunks = line.split("\t");
				Row row = new Row();
				for (int i = 0; i < propertyNames.length; i++) {
					Class<?> propertyType = columns.get(propertyNames[i]).getType();
					if (propertyType == null  || propertyType.equals(String.class)) {
						row.set(propertyNames[i], chunks[i]);
					} else {
						if (propertyType.equals(File.class)) {
							row.set(propertyNames[i], new File(chunks[i]));
						}
						if (propertyType.equals(Integer.class)) {
							row.set(propertyNames[i], Integer.parseInt(chunks[i]));
						}
						if (propertyType.equals(Double.class)) {
							row.set(propertyNames[i], Double.parseDouble(chunks[i]));
						}
						if (propertyType.equals(Boolean.class)) {
							row.set(propertyNames[i], Boolean.parseBoolean(chunks[i]));
						}
					}
				}
				add(row);

				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	/**
	 * 
	 * @return the number of unique rows
	 */
	public int size() {
		return rows.size();
	}
	
	public ArrayList<String> intersectKeyValues(Table other) {
		ArrayList<String> ourKeys = getKeyValues(); 
		
		/* output should be on the order of the size of our rows */
		ArrayList<String> out = new ArrayList<String>(ourKeys.size());
		
		for (String ourKey: ourKeys) {
			if (other.get(ourKey) != null) out.add(ourKey);
		}
		
		return out;
	}
	
	
	/**
	 * Returns a Table that is the combination of this one and another
	 * 
	 * Uses our key as the key column
	 * 
	 * @param other
	 * @return
	 */
	public Table combine(Table other) {
		Table out = new Table(this.getName() + " + " + other.getName(), key);
		ArrayList<Row> rows;
		rows = getRows();
		for (Row row: rows) out.add(row);
		rows = other.getRows();
		for (Row row: rows) out.add(row);
		
		return out;
		
	}

	
}
