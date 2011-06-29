from database import url
from sqlalchemy import create_engine
from razi.tests.common import functions

class TestFunctions(functions.TestFunctions):

    def setup(self):
        self.engine = create_engine(url, echo=True)


    
