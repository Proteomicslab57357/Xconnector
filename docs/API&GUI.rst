API & GUI Implementation
============================

* The API function connects databases in Xconnector is made to be programmatically efficient. Using python generators implementation, only one query is called from the database each time by the API. This will reduce the memory used by Xconnector, as well as overcome the errors that could occur during the slow internet connection.

* After the API sends the output to the GUI. Xconnector utilises multithreading to allow efficient execution for the GUI, which allow multitasking and converting data between the GUI and the API.