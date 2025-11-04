# test_queue.py
import pytest
from queue import Queue

class TestQueue:
    
    def setup_method(self):
        """Setup before each test method"""
        self.queue = Queue()
    
    def teardown_method(self):
        """Cleanup after each test method"""
        self.queue = None
    
    # Enqueue functionality
    def test_enqueue_basic(self):
        
        self.queue.enqueue(1)
        self.queue.enqueue(2)
        assert self.queue.size() == 2
    
    # Dequeue functionality
    def test_dequeue_basic(self):
        
        self.queue.enqueue(1)
        self.queue.enqueue(2)
        result = self.queue.dequeue()
        assert result == 1
        assert self.queue.size() == 1
    
    # Size method
    def test_size(self):
        
        assert self.queue.size() == 0
        self.queue.enqueue("test")
        assert self.queue.size() == 1
        self.queue.enqueue("another")
        assert self.queue.size() == 2
        self.queue.dequeue()
        assert self.queue.size() == 1
    
    # is_empty method
    def test_is_empty(self):
        
        assert self.queue.is_empty() == True
        self.queue.enqueue(1)
        assert self.queue.is_empty() == False
        self.queue.dequeue()
        assert self.queue.is_empty() == True
    
    # Edge case 1: dequeue from empty queue
    def test_dequeue_empty(self):
        """Edge case 1: Test dequeue on empty queue"""
        with pytest.raises(IndexError, match="Cannot dequeue from an empty queue"):
            self.queue.dequeue()
    

    # Edge case 2: different data types
    def test_different_data_types(self):
        """Edge case 2: Test queue with different data types"""
        test_data = [1, "string", 3.14, [1, 2, 3], {"key": "value"}]
        
        for item in test_data:
            self.queue.enqueue(item)
        
        for expected in test_data:
            assert self.queue.dequeue() == expected
    
    #  Edge case 3: very large queue
    def test_large_queue(self):
        """Edge case 3: Test with lots of elements like large queue"""
        for i in range(1000):
            self.queue.enqueue(i)
        
        assert self.queue.size() == 1000
        
        for i in range(1000):
            assert self.queue.dequeue() == i
        
        assert self.queue.is_empty()
    
    # Edge case 4:  Single element
    def test_single_element(self):
        """Edge case 4: Test queue with single element"""
        self.queue.enqueue("single")
        assert self.queue.size() == 1
        assert self.queue.is_empty() == False
        result = self.queue.dequeue()
        assert result == "single"
        assert self.queue.is_empty() == True
