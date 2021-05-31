package UI;

public interface Listener {

    /*
     * When a view's state changes (e.g. a new index of a list is selected)
     * then it tells its Adapter it has changed and passes a copy of itself.
     * The Adapter can then intelligently inquire of the View what was changed
     * and then react accordingly.
     */
    public void viewStateChanged(View updatedView);

}