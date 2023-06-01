package CLPSO_MVNS;

import java.util.AbstractSequentialList;
import java.util.Collection;
import java.util.ConcurrentModificationException;
import java.util.List;
import java.util.ListIterator;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.function.Consumer;

public class MyLinkedList<E> extends AbstractSequentialList<E> implements List<E>,Cloneable{

	int size = 0;
	Node<E> first;
	Node<E> last;
	
	public MyLinkedList(){
	}
	
	public MyLinkedList(Collection<? extends E> c) {
		this();
        addAll(c);
	}
	
	public boolean addAll(Collection<? extends E> c) {
        return addAll(size, c);
    }

    public boolean addAll(int index, Collection<? extends E> c) {
        checkPositionIndex(index);

        Object[] a = c.toArray();
        int numNew = a.length;
        if (numNew == 0)
            return false;

        Node<E> pred, succ;
        if (index == size) {
            succ = null;
            pred = last;
        } else {
            succ = node(index);
            pred = succ.prev;
        }

        for (Object o : a) {
            @SuppressWarnings("unchecked") E e = (E) o;
            Node<E> newNode = new Node<>(pred, e, null);
            if (pred == null)
                first = newNode;
            else
                pred.next = newNode;
            pred = newNode;
        }

        if (succ == null) {
            last = pred;
        } else {
            pred.next = succ;
            succ.prev = pred;
        }

        size += numNew;
        modCount++;
        return true;
    }

	static class Node<E> {
		E item;
		Node<E> next;
		Node<E> prev;

		Node(Node<E> prev, E element, Node<E> next) {
			this.item = element;
			this.next = next;
			this.prev = prev;
		}

		@Override
		public String toString() {
			String s1 = item.toString();
			String s2 = null;
			String s3 = null;
			if (next == null) {
				s2 = "null";
			} else {
				s2 = next.item.toString();
			}
			if (prev == null) {
				s3 = "null";
			} else {
				s3 = prev.item.toString();
			}

			String s = "prev:" + s3 + " item:" + s1 + " next:" + s2;
			return s;
		}
	}

	public void twoOpt(int arc1, int arc2) {
		if (arc1 < 0 || arc2 < 0 || arc1 > size || arc2 > size || arc2 - arc1 < 0
				|| (arc2 - arc1 > 0 && arc2 - arc1 < 2)) {
			throw new NoSuchElementException();
		}
		Node<E> arc2_1 = node(arc2 - 1);
		Node<E> arc2_2 = arc2_1.next;
		Node<E> arc1_2 = null;
		arc2_1.next = arc2_1.prev;
		Node<E> temp = arc2_1.next;
		for (int i = 0; i < arc2 - arc1 - 1; i++) {
			Node<E> temp1 = temp.prev;
			temp.prev = temp.next;
			if (i != arc2 - arc1 - 2) {
				temp.next = temp1;
				temp = temp.next;
			} else {
				temp.next = arc2_2;
				arc2_1.prev = temp1;
				arc1_2 = temp;
				if (temp1 != null) {
					temp1.next = arc2_1;
				}
				if (arc2_2 != null) {
					arc2_2.prev = temp;
				}
			}
		}
		if (arc1 == 0) {
			first = arc2_1;
		}
		if (arc2 == size) {
			last = arc1_2;
		}
	}

	public void orOpt(int len, int begin, int insert) {
		if (begin < 0 || len < 1 || begin + len > size || insert < -1 || insert > size - 1 || (insert >= begin - 1 && insert<=begin+len-1)) {
			throw new NoSuchElementException();
		}
		Node<E> beg = node(begin);
		Node<E> end = null;
		if (len == 1) {
			end = beg;
		} else {
			Node<E> temp = beg.next;
			for (int i = 0; i < len - 2; i++) {
				temp = temp.next;
			}
			end = temp;
		}
		Node<E> endPNext=end.next;
		Node<E> ins=null;
		Node<E> insPNext=null;
		if(insert!=-1) {
			ins = node(insert);
			insPNext = ins.next;
			ins.next = beg;
		}else {
			insPNext=first;
			first=beg;
		}
		Node<E> begPPrev = beg.prev;
		if(endPNext==null) {
			last=begPPrev;
		}
		beg.prev = ins;
		if(endPNext!=null) {
			endPNext.prev = begPPrev;
		}
		if (begPPrev != null) {
			begPPrev.next = endPNext;
		} else {
			first = endPNext;
		}
		end.next = insPNext;
		if (insPNext != null) {
			insPNext.prev = end;
		} else {
			last = end;
		}
	}

	public void addFirst(E e) {
		linkFirst(e);
	}

	private void linkFirst(E e) {
		final Node<E> f = first;
		final Node<E> newNode = new Node<>(null, e, f);
		first = newNode;
		if (f == null)
			last = newNode;
		else
			f.prev = newNode;
		size++;
		modCount++;
	}

	public boolean add(E e) {
		linkLast(e);
		return true;
	}

	void linkLast(E e) {
		final Node<E> l = last;
		final Node<E> newNode = new Node<>(l, e, null);
		last = newNode;
		if (l == null)
			first = newNode;
		else
			l.next = newNode;
		size++;
		modCount++;
	}

	void linkBefore(E e, Node<E> succ) {
		// assert succ != null;
		final Node<E> pred = succ.prev;
		final Node<E> newNode = new Node<>(pred, e, succ);
		succ.prev = newNode;
		if (pred == null)
			first = newNode;
		else
			pred.next = newNode;
		size++;
		modCount++;
	}

	public E remove(int index) {
		checkElementIndex(index);
		return unlink(node(index));
	}

	private void checkElementIndex(int index) {
		if (!isElementIndex(index))
			throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
	}

	private boolean isElementIndex(int index) {
		return index >= 0 && index < size;
	}

	private String outOfBoundsMsg(int index) {
		return "Index: " + index + ", Size: " + size;
	}

	Node<E> node(int index) {
		// assert isElementIndex(index);

		if (index < (size >> 1)) {
			Node<E> x = first;
			for (int i = 0; i < index; i++)
				x = x.next;
			return x;
		} else {
			Node<E> x = last;
			for (int i = size - 1; i > index; i--)
				x = x.prev;
			return x;
		}
	}

	E unlink(Node<E> x) {
		// assert x != null;
		final E element = x.item;
		final Node<E> next = x.next;
		final Node<E> prev = x.prev;

		if (prev == null) {
			first = next;
		} else {
			prev.next = next;
			x.prev = null;
		}

		if (next == null) {
			last = prev;
		} else {
			next.prev = prev;
			x.next = null;
		}

		x.item = null;
		size--;
		modCount++;
		return element;
	}

	@Override
	public ListIterator<E> listIterator(int index) {
		checkPositionIndex(index);
		return new ListItr(index);
	}

	private void checkPositionIndex(int index) {
		if (!isPositionIndex(index))
			throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
	}

	private boolean isPositionIndex(int index) {
		return index >= 0 && index <= size;
	}

	private class ListItr implements ListIterator<E> {
		private Node<E> lastReturned;
		private Node<E> next;
		private int nextIndex;
		private int expectedModCount = modCount;

		ListItr(int index) {
			// assert isPositionIndex(index);
			next = (index == size) ? null : node(index);
			nextIndex = index;
		}

		public boolean hasNext() {
			return nextIndex < size;
		}

		public E next() {
			checkForComodification();
			if (!hasNext())
				throw new NoSuchElementException();

			lastReturned = next;
			next = next.next;
			nextIndex++;
			return lastReturned.item;
		}

		public boolean hasPrevious() {
			return nextIndex > 0;
		}

		public E previous() {
			checkForComodification();
			if (!hasPrevious())
				throw new NoSuchElementException();

			lastReturned = next = (next == null) ? last : next.prev;
			nextIndex--;
			return lastReturned.item;
		}

		public int nextIndex() {
			return nextIndex;
		}

		public int previousIndex() {
			return nextIndex - 1;
		}

		public void remove() {
			checkForComodification();
			if (lastReturned == null)
				throw new IllegalStateException();

			Node<E> lastNext = lastReturned.next;
			unlink(lastReturned);
			if (next == lastReturned)
				next = lastNext;
			else
				nextIndex--;
			lastReturned = null;
			expectedModCount++;
		}

		public void set(E e) {
			if (lastReturned == null)
				throw new IllegalStateException();
			checkForComodification();
			lastReturned.item = e;
		}

		public void add(E e) {
			checkForComodification();
			lastReturned = null;
			if (next == null)
				linkLast(e);
			else
				linkBefore(e, next);
			nextIndex++;
			expectedModCount++;
		}

		public void forEachRemaining(Consumer<? super E> action) {
			Objects.requireNonNull(action);
			while (modCount == expectedModCount && nextIndex < size) {
				action.accept(next.item);
				lastReturned = next;
				next = next.next;
				nextIndex++;
			}
			checkForComodification();
		}

		final void checkForComodification() {
			if (modCount != expectedModCount)
				throw new ConcurrentModificationException();
		}
	}

	@Override
	public int size() {
		return size;
	}
	
	@SuppressWarnings("unchecked")
    private MyLinkedList<E> superClone() {
        try {
            return (MyLinkedList<E>) super.clone();
        } catch (CloneNotSupportedException e) {
            throw new InternalError(e);
        }
    }
	
	public MyLinkedList<E> clone() {
        MyLinkedList<E> clone = superClone();

        // Put clone into "virgin" state
        clone.first = clone.last = null;
        clone.size = 0;
        clone.modCount = 0;

        // Initialize clone with our elements
        for (Node<E> x = first; x != null; x = x.next)
            clone.add(x.item);

        return clone;
    }

}
