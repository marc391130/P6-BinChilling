from io import TextIOWrapper
from time import time
import logging
logging.basicConfig()



class Logger:
    def __init__(self, console_log: bool = True, loggingPath: TextIOWrapper = None) -> None:
        self.__console_log_realtime__ = console_log
        self.__start_time__ = time() 
        self.__Current_Section__ = (1, 0)
        self.loggingPath = loggingPath
    
    #adds new header section to the log
    def AddSection(self, header: str = None) -> None:
        if header is not None:
            self.Log(f"____________________|{header}|____________________")
        
        self.__Current_Section__ = (self.__Current_Section__[0]+1, 0)
        self.StatusLog(f"Runtime: {time() - self.__start_time__}s")
    
    def AddSubsection(self) -> None:
        self.__Current_Section__ = (self.__Current_Section__[0], self.__Current_Section__[1] + 1)
    
    #adds text appended with status including the status
    def StatusLog(self, log: str) -> None:
        self.AddSubsection()
        status = f"Section {self.__Current_Section__}: "
        self.Log(status + log)
    
    #adds text to the logger
    def Log(self, log: str) -> None:
        print(log, file=self.loggingPath)
    
    #logs an error
    def LogError(self, exception: Exception) -> None:
        self.log(str(exception))
    
    def AddSpace(self) -> None:
        self.log('\n')
        
    