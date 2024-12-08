import multiprocessing
import threading
import time

def thread_task1(instance_id, arg1, arg2):
    print(f"Thread 1 starting for Instance {instance_id} with args ({arg1}, {arg2}) (PID: {multiprocessing.current_process().pid})")
    time.sleep(2)  # Simulate work
    print(f"Thread 1 finished for Instance {instance_id} with args ({arg1}, {arg2}) (PID: {multiprocessing.current_process().pid})")

def thread_task2(instance_id, arg1, arg2):
    print(f"Thread 2 starting for Instance {instance_id} with args ({arg1}, {arg2}) (PID: {multiprocessing.current_process().pid})")
    time.sleep(3)  # Simulate work
    print(f"Thread 2 finished for Instance {instance_id} with args ({arg1}, {arg2}) (PID: {multiprocessing.current_process().pid})")

def program1(args):
    # Unpack the arguments
    instance_id, arg1, arg2 = args
    print(f"Program1 Instance {instance_id} starting with args ({arg1}, {arg2}) (PID: {multiprocessing.current_process().pid})")
    
    # Create and start threads
    thread1 = threading.Thread(target=thread_task1, args=(instance_id, arg1, arg2))
    thread2 = threading.Thread(target=thread_task2, args=(instance_id, arg1, arg2))
    
    thread1.start()
    thread2.start()
    
    # Wait for threads to finish
    thread1.join()
    thread2.join()
    print(f"Program1 Instance {instance_id} finished with args ({arg1}, {arg2}) (PID: {multiprocessing.current_process().pid})")

def main():
    total_instances = 20
    num_processes = 12  # Number of processes to run concurrently

    # Create a list of tasks with instance_id and the additional arguments
    tasks = [(i, f"arg1_{i}", f"arg2_{i}") for i in range(total_instances)]
    
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Map the program1 function to the tasks
        pool.map(program1, tasks)

if __name__ == "__main__":
    main()