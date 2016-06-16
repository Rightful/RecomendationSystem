
public class Movie {
    
    int index, year;
    String title;
    
	public Movie(int _index, int _year, String _title) {
        this.index = _index;
        this.year  = _year;
        this.title = _title;
    }
    
    public int getIndex() {
        return index;
    }
    
    public int getYear() {
        return year;
    }
    
    public String getTitle() {
        return title;
    }
    
}

