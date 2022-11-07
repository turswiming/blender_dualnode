# ***** BEGIN GPL LICENSE BLOCK *****
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ***** END GPL LICENSE BLOCK *****

if(WIN32)
  set(XML2_EXTRA_ARGS 
    -DLIBXML2_WITH_ZLIB=OFF
    -DLIBXML2_WITH_LZMA=OFF
    -DLIBXML2_WITH_PYTHON=OFF
    -DLIBXML2_WITH_ICONV=OFF
    -DLIBXML2_WITH_TESTS=OFF
    -DLIBXML2_WITH_PROGRAMS=OFF
    -DBUILD_SHARED_LIBS=OFF
  )
  ExternalProject_Add(external_xml2
    URL file://${PACKAGE_DIR}/${XML2_FILE}
    DOWNLOAD_DIR ${DOWNLOAD_DIR}
    URL_HASH ${XML2_HASH_TYPE}=${XML2_HASH}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${LIBDIR}/xml2 ${DEFAULT_CMAKE_FLAGS} ${XML2_EXTRA_ARGS}
    PREFIX ${BUILD_DIR}/xml2
    INSTALL_DIR ${LIBDIR}/xml2
  )
else()
  ExternalProject_Add(external_xml2
    URL file://${PACKAGE_DIR}/${XML2_FILE}
    DOWNLOAD_DIR ${DOWNLOAD_DIR}
    URL_HASH ${XML2_HASH_TYPE}=${XML2_HASH}
    PREFIX ${BUILD_DIR}/xml2
    CONFIGURE_COMMAND ${CONFIGURE_ENV} && cd ${BUILD_DIR}/xml2/src/external_xml2/ && ${CONFIGURE_COMMAND}
      --prefix=${LIBDIR}/xml2
      --disable-shared
      --enable-static
      --with-pic
      --with-python=no
      --with-lzma=no
      --with-zlib=no
      --with-iconv=no
    BUILD_COMMAND ${CONFIGURE_ENV} && cd ${BUILD_DIR}/xml2/src/external_xml2/ && make -j${MAKE_THREADS}
    INSTALL_COMMAND ${CONFIGURE_ENV} && cd ${BUILD_DIR}/xml2/src/external_xml2/ && make install
    INSTALL_DIR ${LIBDIR}/xml2
  )
endif()

if(WIN32 AND BUILD_MODE STREQUAL Release)
  ExternalProject_Add_Step(external_xml2 after_install
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${LIBDIR}/xml2/include  ${HARVEST_TARGET}/xml2/include
        COMMAND ${CMAKE_COMMAND} -E copy ${LIBDIR}/xml2/lib/libxml2s.lib  ${HARVEST_TARGET}/xml2/lib/libxml2s.lib
    DEPENDEES install
  )
endif()
