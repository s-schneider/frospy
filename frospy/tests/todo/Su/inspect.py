from frospy.preprocessing.inspect import inspect
import sys

#event = sys.argv[1] # Event name
#event='060994A'
event='100494B'
inspect('300hrs/raw/%s'%event, outputformat='pickle', raw_only=False, start=3)
