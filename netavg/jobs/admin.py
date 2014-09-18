from django.contrib import admin
from .models import Job, Result


class JobAdmin(admin.ModelAdmin):
    pass


class ResultAdmin(admin.ModelAdmin):
    pass


''' Register Admin layouts into django'''
admin.site.register(Job, JobAdmin)
admin.site.register(Result, ResultAdmin)
