//*** BEGIN XUE MINGQIANG **** ///
import java.util.*;

class Node {
    private int rank;
    private ArrayList<Node> outgoing;
    private Node preImage;
    private Node postImage;
    private Node child;
    private ArrayList<Node> visited;
    private boolean inCycle = false;


    public Node(){
        outgoing = new ArrayList<Node>();
        visited = new ArrayList<Node>();
    }

    public boolean hasVisited(Node next){
        return visited.contains(next);
    }


    public boolean isInCycle(){
        return inCycle;
    }

    public void setInCycle(boolean inCycle){
        this.inCycle = inCycle;
    }
    
    public void addToVisited(Node next){
        visited.add(next);
    }

    public void resetVisited(){
        visited = new ArrayList<Node>();
    }

    public Node getPreImage(){
        return preImage;
    }

    public Node getPostImage(){
        return postImage;
    }

    public void setPreImage(Node preImage){
        this.preImage = preImage;
    }

    public void setPostImage(Node postImage){
        this.postImage = postImage;
    }

    public void addOutgoing(Node link){
        outgoing.add(link);
    }

    public ArrayList<Node> getAllOutgoing(){
        return outgoing;
    }

    public int getRank(){
        return rank;
    }

    public void setRank(int rank){
        this.rank = rank;
    }

    public void setChild(Node child){
        this.child = child;
    }

    public Node getChild(){
        return child;
    }
}

public class Randomization {
    public static Random rand = new Random(10); //seeded by 10

    public static Node selectNext(Node t_i, ArrayList<Node> nodeList){
        Node _t_j = null;
        boolean isFound = false;
        Node lastChoice = null;
        boolean foundLastChoice = false;
        ArrayList<Node> outgoing = t_i.getAllOutgoing();
        int randOffset = rand.nextInt(outgoing.size());
        int idx = randOffset;
        for(int i=0; i<outgoing.size(); i++) {
            _t_j = outgoing.get(idx % outgoing.size());
            if(t_i.hasVisited(_t_j)){
                idx++;
                continue;
            }
            if(!foundLastChoice && _t_j.isInCycle()){
                lastChoice = _t_j;
                foundLastChoice = true;
            }
            if(!_t_j.isInCycle()){
                isFound = true;
                break;
            }
            idx++;
        }
        if(!isFound)
            return lastChoice;
        return _t_j;
    }

    public static int [] getRandAssignment(ArrayList<Node> nodeList){
        ArrayList<Node> L = new ArrayList<Node>(); // the list of unprocessed nodes

        // initialize unprocessed
        for(int i=0; i<nodeList.size(); i++)
            L.add(nodeList.get(i));
        while(!L.isEmpty()){
            Node t_i = L.get(rand.nextInt(L.size()));
            Node _t_i = t_i.getPostImage();
            Node _t_j = selectNext(t_i, nodeList);
            if(_t_j == null){
                System.out.println("Error: _t_j == null");
                System.exit(1);
            }
            Node temp = _t_j;
            ArrayList<Node> walkTrack = new ArrayList<Node>(); // track all nodes in the walk

            t_i.addToVisited(_t_j);
            _t_j.setInCycle(true);
            walkTrack.add(t_i);
            while(_t_j!=_t_i){
                Node t_x = _t_j.getPreImage();
                Node _t_y = selectNext(t_x, nodeList);
                if(_t_y == null){
                    System.out.println("Error: _t_j == null");
                    System.exit(1);
                }
                t_x.addToVisited(_t_y);
                walkTrack.add(t_x);
                _t_y.setInCycle(true);
                _t_j.setChild(_t_y);
                _t_j = _t_y;
            }
            _t_i.setChild(temp);
            Node next = _t_i;
            ArrayList<Node> preTrack = new ArrayList<Node>();
            ArrayList<Node> postTrack = new ArrayList<Node>();
            while(true){
                Node pre = next.getPreImage();
                Node post = next.getChild();
                preTrack.add(pre);
                postTrack.add(post);
                next = next.getChild();
                if(next == _t_i)
                    break;
            }
            L.removeAll(preTrack);

            for(int i=0; i<walkTrack.size(); i++){
                walkTrack.get(i).resetVisited();
                walkTrack.get(i).getPostImage().setInCycle(false);
            }

            for(int i=0; i<preTrack.size(); i++){
                Node pre = preTrack.get(i);
                Node post = postTrack.get(i);
                int preRank = pre.getRank();
                int postRank = post.getRank();
                nodeList.get(preRank).setPostImage(post);
                nodeList.get(postRank).setPreImage(pre);
            }
        }

        int [] A = new int[nodeList.size()];
        for(int i=0; i<nodeList.size(); i++){
            Node pre = nodeList.get(i);
            Node post = pre.getPostImage();
            A[pre.getRank()] = post.getRank();
        }
        return A;
    }

    public static int [] run(int [][] final_assignment, int offset, int size, int k){
        ArrayList<Node> nodeList = new ArrayList<Node>();

        //initialize the nodes
        for(int i=0; i<size; i++){
            Node node = new Node();
            node.setRank(i); // the id of a node is the index of a record in dataset
            node.setPreImage(node);
            node.setPostImage(node);
            nodeList.add(node);
        }

        // update the matchings
        for(int i=0; i<k; i++){
            //System.out.println("round " + i +":");
            for (int j=0; j<size; j++){
                int match = final_assignment[j][i];
                nodeList.get(j).getAllOutgoing().add(nodeList.get(match));
            }
        }
        System.out.println("Random assignment created");
        int [] A = getRandAssignment(nodeList);
        /*
        for(int i=0; i<A.length; i++){
            System.out.println(i + " " + A[i]);
        }
        */
        return A;
    }

}
//*** END XUE MINGQIANG **** ///

