package Database;

/**
 * This does not contain a column of data.  More, it is the
 * description of the column.
 * <p>
 * Copyright 2013, Brian Risk
 *
 * @author Brian Risk
 */
public class Column {

    private String name;
    private Class<?> type;
    private Object defaultValue;

    public Column(String name, Class<?> type, Object defaultValue) {
        super();
        this.name = name;
        this.type = type;
        this.defaultValue = defaultValue;
    }

    public Column(String name, Class<?> type) {
        this(name, type, null);
    }

    public String getName() {
        return name;
    }

    public Class<?> getType() {
        return type;
    }

    public Object getDefaultValue() {
        return defaultValue;
    }

    public boolean equals(Column other) {
        if (!name.equals(other.getName())) return false;
//		if (!type.equals(other.getType())) return false;
        return true;
    }


}
