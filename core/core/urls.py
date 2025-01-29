"""core URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.10/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.urls import include, re_path
from django.contrib import admin
from django.conf import settings
from django.conf.urls.static import static


urlpatterns = [
    re_path(r'^admin/', admin.site.urls),
    re_path('accounts/', include('accounts.urls', namespace='accounts')),
    re_path('organisms/', include('organisms.urls', namespace='organisms')),
    re_path('predictors/', include('predictors.urls', namespace='predictors')),
    re_path('h5addatasets/', include('h5addatasets.urls', namespace='h5addatasets')),
    re_path('rdsdatasets/', include('rdsdatasets.urls', namespace='rdsdatasets')),
    re_path('tenxfeaturebcmatrixdatasets/', include('tenxfeaturebcmatrixdatasets.urls', namespace='tenxfeaturebcmatrixdatasets')),
    re_path('preddatasets/', include('preddatasets.urls', namespace='preddatasets')),
    re_path('scripts/', include('scripts.urls', namespace='scripts')),
    re_path('jobs/', include('jobs.urls', namespace='jobs')),
]

urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
