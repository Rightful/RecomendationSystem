
public class Rating {

    User user;
    Movie movie;
    double rating;
    
	public Rating(User _user, Movie _movie, int _rating) {
        this.user      = _user;
        this.movie     = _movie;
        this.rating    = (double) _rating;
    }
    
    public Rating(User _user, Movie _movie, double _rating) {
        this.user      = _user;
        this.movie     = _movie;
        this.rating    = _rating;
    }

    
    public User getUser() {
        return user;
    }
    
    public Movie getMovie() {
        return movie;
    }
    
    public double getRating() {
        return rating;
    }
    
    public void setRating(double _rating) {
        this.rating = _rating;
    }
}
