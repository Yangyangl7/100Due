import java.util.*;
import java.util.AbstractMap.SimpleEntry;
import java.io.*; 
import java.awt.*;

public class Assignment2{
	class Point {
    	public int x;
    	public int y;

    	Point(int x, int y) {
        	this.x = x;
        	this.y = y;
    	}

    	public int crossProduct(Point p) {
        	return x * p.y - p.x * y;
    	}

    	public Point substract(Point p) {
        	return new Point(x - p.x, y - p.y);
    	}

    	public static int direction(Point pi, Point pj, Point pk) {
        	return pk.substract(pi).crossProduct(pj.substract(pi));
    	}

    	public static boolean onSegment(Point pi, Point pj, Point pk) {
        	if (Math.min(pi.x, pj.x) <= pk.x
                	&& pk.x <= Math.max(pi.x, pj.x)
                	&& Math.min(pi.y, pj.y) <= pk.y
                	&& pk.y <= Math.max(pi.y, pj.y))
            	return true;
        	return  false;
    	}
	}
	class EndPoint extends Point implements Comparable<EndPoint> {
		public Segment segment;
	
		public EndPoint(int x, int y) {
			super(x, y);
			this.segment = null;
		}
	
		public boolean isLeft() {
			return this.segment != null && this.segment.left == this;
		}
	
		@Override
		public int compareTo(EndPoint o) {
			if (x < o.x) {
				return -1;
			} else if (x > o.x) {
				return 1;
			} else if(this.isLeft() && !o.isLeft()) {
				return -1;
			} else if (!this.isLeft() && o.isLeft()) {
				return 1;
			} else if (y < o.y) {
				return -1;
			} else if (y > o.y) {
				return 1;
			}
			return 0;
		}
	}
	
	class Segment implements Comparable<Segment> {
		public EndPoint left;
		public EndPoint right;
		SweepLineStatus T;
	
		Segment(EndPoint p1, EndPoint p2) {
			if (p1.compareTo(p2) > 0) {
				this.left = p2;
				this.right = p1;
			} else {
				this.left = p1;
				this.right = p2;
			}
			this.left.segment = this;
			this.right.segment = this;
			this.T = null;
		}
	
		public boolean intersect(Segment s) {
			EndPoint p1 = this.left;
			EndPoint p2 = this.right;
			EndPoint p3 = s.left;
			EndPoint p4 = s.right;
	
			double d1 = Point.direction(p3,p4,p1);
			double d2 = Point.direction(p3,p4,p2);
			double d3 = Point.direction(p1,p2,p3);
			double d4 = Point.direction(p1,p2,p4);
			if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0))
					&& ((d3 > 0 && d4 < 0)||(d3 < 0 && d4 > 0)))
				return true;
			else if(d1 == 0 && Point.onSegment(p3,p4,p1))
				return true;
			else if(d2 == 0 && Point.onSegment(p3,p4,p2))
				return true;
			else if(d3 == 0 && Point.onSegment(p1,p2,p3))
				return true;
			else if(d4 ==0 && Point.onSegment(p1,p2,p4))
				return true;
			return false;
		}
	
		@Override
		public int compareTo(Segment o) {
			if (intersect(o)) {
	
			} else {
	
			}
		}
	}
	// Heapsort for sorting segments
	class Heapsort<AnyType extends Comparable<? super AnyType>> {
		public void sort(AnyType arr[]) {
			int n = arr.length;
			// Build heap (rearrange array)
			for (int i = n / 2 - 1; i >= 0; i--) {
				heapify(arr, n, i);
			}
			// One by one extract an element from heap
			for (int i = n - 1; i > 0; i--) {
				// Move current root to end
				AnyType temp = arr[0];
				arr[0] = arr[i];
				arr[i] = temp;
				// call max heapify on the reduced heap
				heapify(arr, i, 0);
			}
		}
	
		// Heapify a subtree
		void heapify(AnyType[] arr, int n, int i) {
			int largest = i; // Initialize largest as root
			int l = 2 * i + 1; // left = 2*i + 1
			int r = 2 * i + 2; // right = 2*i + 2
			// If left child is larger than root
			if (l < n && arr[l].compareTo(arr[largest]) > 0) {
				largest = l;
			}
			// If right child is larger than largest so far
			if (r < n && arr[r].compareTo(arr[largest]) > 0) {
				largest = r;
			}
			// If largest is not root
			if (largest != i) {
				AnyType swap = arr[i];
				arr[i] = arr[largest];
				arr[largest] = swap;
				// Recursively heapify the affected sub-tree
				heapify(arr, n, largest);
			}
		}
	}	
	// Use self-adjust BST for maintaining status
	class AvlTree<AnyType extends Comparable<? super AnyType>> {
		/**
		 * Construct the tree.
		 */
		public AvlTree() {
			root = null;
		}
	
		/**
		 * Insert into the tree; duplicates are ignored.
		 *
		 * @param x the item to insert.
		 */
		public void insert(AnyType x) {
			root = insert(x, root);
		}
	
		/**
		 * Remove from the tree. Nothing is done if x is not found.
		 *
		 * @param x the item to remove.
		 */
		public void remove(AnyType x) {
			root = remove(x, root);
		}
	
		/**
		 * Internal method to remove from a subtree.
		 *
		 * @param x the item to remove.
		 * @param t the node that roots the subtree.
		 * @return the new root of the subtree.
		 */
		private AvlNode<AnyType> remove(AnyType x, AvlNode<AnyType> t) {
			if (t == null)
				return t;   // Item not found; do nothing
	
			int compareResult = x.compareTo(t.element);
	
			if (compareResult < 0)
				t.left = remove(x, t.left);
			else if (compareResult > 0)
				t.right = remove(x, t.right);
			else if (t.left != null && t.right != null) // Two children
			{
				t.element = findMin(t.right).element;
				t.right = remove(t.element, t.right);
			} else
				t = (t.left != null) ? t.left : t.right;
			return balance(t);
		}
	
		/**
		 * Find the smallest item in the tree.
		 *
		 * @return smallest item or null if empty.
		 * @throws Exception
		 */
		public AnyType findMin() throws Exception {
			if (isEmpty())
				throw new Exception();
			return findMin(root).element;
		}
	
		/**
		 * Find the largest item in the tree.
		 *
		 * @return the largest item of null if empty.
		 * @throws Exception
		 */
		public AnyType findMax() throws Exception {
			if (isEmpty())
				throw new Exception();
			return findMax(root).element;
		}
	
		/**
		 * Find an item in the tree.
		 *
		 * @param x the item to search for.
		 * @return true if x is found.
		 */
		public boolean contains(AnyType x) {
			return contains(x, root);
		}
	
		/**
		 * Make the tree logically empty.
		 */
		public void makeEmpty() {
			root = null;
		}
	
		/**
		 * Test if the tree is logically empty.
		 *
		 * @return true if empty, false otherwise.
		 */
		public boolean isEmpty() {
			return root == null;
		}
	
		/**
		 * Print the tree contents in sorted order.
		 */
		public void printTree() {
			if (isEmpty())
				System.out.println("Empty tree");
			else
				printTree(root);
		}
	
		private static final int ALLOWED_IMBALANCE = 1;
	
		// Assume t is either balanced or within one of being balanced
		private AvlNode<AnyType> balance(AvlNode<AnyType> t) {
			if (t == null)
				return t;
			if (height(t.left) - height(t.right) > ALLOWED_IMBALANCE)
				if (height(t.left.left) >= height(t.left.right))
					t = rotateWithLeftChild(t);
				else
					t = doubleWithLeftChild(t);
			else if (height(t.right) - height(t.left) > ALLOWED_IMBALANCE)
				if (height(t.right.right) >= height(t.right.left))
					t = rotateWithRightChild(t);
				else
					t = doubleWithRightChild(t);
			t.height = Math.max(height(t.left), height(t.right)) + 1;
			return t;
		}
	
		public void checkBalance() {
			checkBalance(root);
		}
	
		private int checkBalance(AvlNode<AnyType> t) {
			if (t == null)
				return -1;
	
			if (t != null) {
				int hl = checkBalance(t.left);
				int hr = checkBalance(t.right);
				if (Math.abs(height(t.left) - height(t.right)) > 1 ||
						height(t.left) != hl || height(t.right) != hr)
					System.out.println("OOPS!!");
			}
			return height(t);
		}
	/**
     * Internal method to insert into a subtree.
     *
     * @param x the item to insert.
     * @param t the node that roots the subtree.
     * @return the new root of the subtree.
     */
    	private AvlNode<AnyType> insert(AnyType x, AvlNode<AnyType> t) {
        	if (t == null)
            	return new AvlNode<AnyType>(x, null, null);

        	int compareResult = x.compareTo(t.element);

        	if (compareResult < 0)
            	t.left = insert(x, t.left);
        	else if (compareResult > 0)
            	t.right = insert(x, t.right);
        	else
            	;  // Duplicate; do nothing
        	return balance(t);
    	}
    /**
     * Internal method to find the smallest item in a subtree.
     *
     * @param t the node that roots the tree.
     * @return node containing the smallest item.
     */
    	private AvlNode<AnyType> findMin(AvlNode<AnyType> t) {
        	if (t == null)
            	return t;

        	while (t.left != null)
           		t = t.left;
        	return t;
    	}
    /**
     * Internal method to find the largest item in a subtree.
     *
     * @param t the node that roots the tree.
     * @return node containing the largest item.
     */
    	private AvlNode<AnyType> findMax(AvlNode<AnyType> t) {
        	if (t == null)
            	return t;

        	while (t.right != null)
            	t = t.right;
        	return t;
    	}
    /**
     * Internal method to find an item in a subtree.
     *
     * @param x is item to search for.
     * @param t the node that roots the tree.
     * @return true if x is found in subtree.
     */
    	private boolean contains(AnyType x, AvlNode<AnyType> t) {
        	while (t != null) {
            	int compareResult = x.compareTo(t.element);

            	if (compareResult < 0)
                	t = t.left;
            	else if (compareResult > 0)
                	t = t.right;
            	else
                	return true;    // Match
        	}
        	return false;   // No match
    	}
    /**
     * Internal method to print a subtree in sorted order.
     *
     * @param t the node that roots the tree.
     */
    	private void printTree(AvlNode<AnyType> t) {
        	if (t != null) {
            	printTree(t.left);
            	System.out.println(t.element);
            	printTree(t.right);
        	}
    	}
    /**
     * Return the height of node t, or -1, if null.
     */
    	private int height(AvlNode<AnyType> t) {
        	return t == null ? -1 : t.height;
    	}
    /**
     * Rotate binary tree node with left child.
     * For AVL trees, this is a single rotation for case 1.
     * Update heights, then return new root.
     */
    	private AvlNode<AnyType> rotateWithLeftChild(AvlNode<AnyType> k2) {
        	AvlNode<AnyType> k1 = k2.left;
       		k2.left = k1.right;
        	k1.right = k2;
        	k2.height = Math.max(height(k2.left), height(k2.right)) + 1;
        	k1.height = Math.max(height(k1.left), k2.height) + 1;
        	return k1;
    	}
    /**
     * Rotate binary tree node with right child.
     * For AVL trees, this is a single rotation for case 4.
     * Update heights, then return new root.
     */
    	private AvlNode<AnyType> rotateWithRightChild(AvlNode<AnyType> k1) {
        	AvlNode<AnyType> k2 = k1.right;
        	k1.right = k2.left;
        	k2.left = k1;
        	k1.height = Math.max(height(k1.left), height(k1.right)) + 1;
        	k2.height = Math.max(height(k2.right), k1.height) + 1;
        	return k2;
    	}

    /**
     * Double rotate binary tree node: first left child
     * with its right child; then node k3 with new left child.
     * For AVL trees, this is a double rotation for case 2.
     * Update heights, then return new root.
     */
    	private AvlNode<AnyType> doubleWithLeftChild(AvlNode<AnyType> k3) {
        	k3.left = rotateWithRightChild(k3.left);
        	return rotateWithLeftChild(k3);
    	}

    /**
     * Double rotate binary tree node: first right child
     * with its left child; then node k1 with new right child.
     * For AVL trees, this is a double rotation for case 3.
     * Update heights, then return new root.
     */
    	private AvlNode<AnyType> doubleWithRightChild(AvlNode<AnyType> k1) {
        	k1.right = rotateWithLeftChild(k1.right);
        	return rotateWithRightChild(k1);
    	}

    	private class AvlNode<AnyType> {
        	// Constructors
        	AvlNode(AnyType theElement) {
            	this(theElement, null, null);
        	}

        	AvlNode(AnyType theElement, AvlNode<AnyType> lt, AvlNode<AnyType> rt) {
            	element = theElement;
            	left = lt;
           		right = rt;
            	height = 0;
        	}
        	AnyType element;      // The data in the node
        	AvlNode<AnyType> left;         // Left child
        	AvlNode<AnyType> right;        // Right child
        	int height;       // Height
    	}
    /**
     * The tree root.
     */
    	private AvlNode<AnyType> root;
	}
	// Radix sort for segments sorting
	// A utility function to get maximum value in arr[]
    static int getMax(int arr[], int n)
    {
        int mx = arr[0];
        for (int i = 1; i < n; i++)
            if (arr[i] > mx)
                mx = arr[i];
        return mx;
    } 
    // A function to do counting sort of arr[] according to
    // the digit represented by exp.
    static void countSort(int arr[], int n, int exp)
    {
        int output[] = new int[n]; // output array
        int i;
        int count[] = new int[10];
        Arrays.fill(count, 0);
 
        // Store count of occurrences in count[]
        for (i = 0; i < n; i++)
            count[(arr[i] / exp) % 10]++;
 
        // Change count[i] so that count[i] now contains
        // actual position of this digit in output[]
        for (i = 1; i < 10; i++)
            count[i] += count[i - 1];
 
        // Build the output array
        for (i = n - 1; i >= 0; i--) {
            output[count[(arr[i] / exp) % 10] - 1] = arr[i];
            count[(arr[i] / exp) % 10]--;
        }
 
        // Copy the output array to arr[], so that arr[] now
        // contains sorted numbers according to current digit
        for (i = 0; i < n; i++)
            arr[i] = output[i];
    }
    // Radix Sort
    public class RadixSort{
    	// A utility function to get maximum value in arr[]
        int getMax(int arr[], int n)
        {
            int mx = arr[0];
            for (int i = 1; i < n; i++)
                if (arr[i] > mx)
                    mx = arr[i];
            return mx;
        }   
        // A function to do counting sort of arr[] according to
        // the digit represented by exp.
        void countSort(int arr[], int n, int exp){
            int output[] = new int[n]; // output array
            int i;
            int count[] = new int[10];
            Arrays.fill(count, 0); 
            // Store count of occurrences in count[]
            for (i = 0; i < n; i++)
                count[(arr[i] / exp) % 10]++;     
            // Change count[i] so that count[i] now contains
            // actual position of this digit in output[]
            for (i = 1; i < 10; i++)
                count[i] += count[i - 1];
            // Build the output array
            for (i = n - 1; i >= 0; i--) {
                output[count[(arr[i] / exp) % 10] - 1] = arr[i];
                count[(arr[i] / exp) % 10]--;
            }
     
            // Copy the output array to arr[], so that arr[] now
            // contains sorted numbers according to current digit
            for (i = 0; i < n; i++)
                arr[i] = output[i];
        }
        void radixsort(int arr[], int n){
            // Find the maximum number to know number of digits
            int m = getMax(arr, n);
     
            // Do counting sort for every digit. Note that
            // instead of passing digit number, exp is passed.
            // exp is 10^i where i is current digit number
            for (int exp = 1; m / exp > 0; exp *= 10)
                countSort(arr, n, exp);
        }
    }
    // vEB Tree for status maintainance
    public class vEBtree <T extends Comparable<T>>{
    	class Node {
            public int u;
            public SimpleEntry<T, T> min;
            public SimpleEntry<T, T> max;
            public Node summary;
            public Node[] cluster;
            public Node(int u) {
                this.u = u;
                min = NIL;
                max = NIL;
                initialize(u);
            }
            private void initialize(int u) {
                if (u <= 2) {
                    summary = null;
                    cluster = null;
                } else {
                    int size = higherSquareRoot();
                    summary = new Node(size);
                    cluster = new vEBtree.Node[size];
                    for (int i = 0; i < size; i++) {
                        cluster[i] = new Node(size);
                    }
                }
            }
    // Deal with higher size value
            private int higherSquareRoot() {
                return (int) Math.pow(2, Math.ceil((Math.log10(u) / Math.log10(2)) / 2));
            }
        }
        // Define NIL value to initialize max;
        private SimpleEntry<T, T> NIL;
        private T ONE;
        private T ZERO;
        private Node root;
        /*
        * Construction method
        */
        public vEBtree(int u, T NIL, T ONE, T ZERO) throws Exception {
            this.NIL = new SimpleEntry<T, T>(NIL, NIL);
            this.ONE = ONE;
            this.ZERO = ZERO;
            if (!isPowerOf2(u)) {
                throw new Exception("Tree size must be a power of 2!");
            }
            root = new Node(u);
        }
        public void insert(T value, T priority) {
            insert(root, value, priority);
        }
        public boolean decreaseKey(T value, T priority) {
            return decreaseKey(root, new SimpleEntry<T, T>(value, priority));
        }
    // Extract maximum value
        public SimpleEntry<T, T> extractMax() {
            SimpleEntry<T, T> max = root.max;
            decreaseKey(max.getKey(), max.getValue());
            return max;
        }
    // Check if the node is null
        private void insertEmptyNode(Node v, SimpleEntry<T, T> s) {
            v.min = s;
            v.max = s;
        }
        /**
        * Insert new value to vEBTree
        * @param v Root
        * @param priority value to insert
        */
        private void insert(Node v, T value, T priority) {
            SimpleEntry<T, T> x = new SimpleEntry<T, T>(value, priority);
            if (v.min == null) {
                insertEmptyNode(v, x);
                return;
            }
            increaseKey(v, x);
        }
        private void increaseKey(Node v, SimpleEntry<T, T> x) {
            if (compareTo(x, v.min) < 0) {
                // Exchange x with v.min
                SimpleEntry<T, T> temp = x;
                x = v.min;
                v.min = temp;
            }
            if (v.u > 2) {
                if (v.cluster[(int) high(v, x)].min == null) {
                    insert(v.summary, x.getKey(), high(v, x));
                    insertEmptyNode(v.cluster[(int) high(v, x)], new SimpleEntry<T, T>(x.getKey(), low(v, x)));
                } else {
                    insert(v.cluster[(int) high(v, x)], x.getKey(), low(v, x));
                }
            }


            if (compareTo(x, v.max) > 0) {
                v.max = x;
            }
        }
        /**
        * Compare 2 elements of heap
        *
        * @param min
        * @param max
        * @return
        */
        private int compareTo(SimpleEntry<T, T> min, SimpleEntry<T, T> max) {


            if (min.getValue().equals(max.getValue())) {
                return min.getKey().compareTo(max.getKey());
            }
            return min.getValue().compareTo(max.getValue());
        }
        /**
        * Compare 2 elements of heap equals or not
        *
        * @param x
        * @param min
        * @return
        */
        private boolean equals(SimpleEntry<T, T> x, SimpleEntry<T, T> min) {
            if (x == null || min == null) {
                return false;
            }
            return x.getKey().equals(min.getKey()) && x.getValue().equals(min.getValue());
        }
        private boolean decreaseKey(Node v, SimpleEntry<T, T> x) {
            if (compareTo(v.min, v.max) == 0) {
                v.min = NIL;
                v.max = NIL;
                return false;
            }
            if (v.u == 2) {
                v.min = ZERO.equals(x.getValue()) ? new SimpleEntry<T, T>(x.getKey(), ONE)
                        : new SimpleEntry<T, T>(x.getKey(), ZERO);
                v.max = v.min;
                return false;
            }
            if (!equals(x, v.min)) {
                return false;
            }
            SimpleEntry<T, T> first_cluster = v.summary.min;
            T priority = index(v, first_cluster, v.cluster[(int) first_cluster.getValue()].min);
            v.min = new SimpleEntry<T, T>(x.getValue(), priority);
            decreaseKey(v.cluster[(int) high(v, x)], new SimpleEntry<T, T>(x.getKey(), low(v, x)));
            if (v.cluster[(int) high(v, x)].min == null) {
                decreaseKey(v.summary, new SimpleEntry<T, T>(x.getKey(), high(v, x)));
                if (equals(x, v.max)) {
                    SimpleEntry<T, T> summary_max = v.summary.max;
                    if (summary_max == null) {
                        v.max = v.min;
                    } else {
                        priority = index(v, summary_max, v.cluster[(int) summary_max.getValue()].max);
                        v.max = new SimpleEntry<T, T>(x.getValue(), priority);
                    }
                }
            } else if (equals(x, v.max)) {
                priority = index(v, new SimpleEntry<T, T>(x.getValue(), high(v, x)), v.cluster[(int) high(v, x)].max);
                v.max = new SimpleEntry<T, T>(x.getValue(), priority);
            }
            return true;
        }
        /*
        * Returns the integer value of the first half of the bits of x.
        */
        @SuppressWarnings({ "deprecation", "unchecked" })
    	private T high(Node node, SimpleEntry<T, T> x) {
            return (T) new Integer((int) Math.floor((int) x.getValue() / lowerSquareRoot(node)));
        }
        /**
        * The integer value of the second half of the bits of x.
        */
        @SuppressWarnings({ "deprecation", "unchecked" })
    	private T low(Node node, SimpleEntry<T, T> x) {
            return (T) new Integer((int) x.getValue() % (int) lowerSquareRoot(node));
        }
        /**
        * The value of the least significant bits of x.
        */
        private double lowerSquareRoot(Node node) {
            return Math.pow(2, Math.floor((Math.log10(node.u) / Math.log10(2)) / 2.0));
        }
        /**
        * The index in the tree of the given value.
        */
        @SuppressWarnings({ "deprecation", "unchecked" })
    	private T index(Node node, SimpleEntry<T, T> first_cluster, SimpleEntry<T, T> min) {
            return (T) new Integer((int) ((int) first_cluster.getValue() * lowerSquareRoot(node) + (int) min.getValue()));
        }
        /**
        * Returns true if x is a power of 2, false otherwise.
        */
        private boolean isPowerOf2(int x) {
            if (0 == x) {
                return false;
            }
            while (x % 2 == 0) {
                x = x / 2;
            }
            if (x > 1) {
                return false;
            }
            return true;
        }
    }
    
    //Sweeping line algorithm 1 (Heapsort and BST)
    
    //Sweeping line algorithm 2 (Radixsort and vEB Tree)
    
    // Use u to determine algorithm 1 or 2
    
    // The main function to compute the running time    
	public void main(String[] args) {


	}
}
