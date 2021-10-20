from django.utils.deprecation import MiddlewareMixin

class HttpResponseCustomHeader(MiddlewareMixin):
    def process_response(self, request, response):
        if not response.has_header("Version"):
            response['Access-Control-Allow-Origin'] = "http://sysu-software.com"
            response["Version"] = "2"
            response["Version_old"] = "1"
        return response
