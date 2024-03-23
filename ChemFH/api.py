from ninja import NinjaAPI, Schema
from service.api import router as service_router

api = NinjaAPI(title='ChemFH')

api.add_router("/", service_router, tags=["service"])
