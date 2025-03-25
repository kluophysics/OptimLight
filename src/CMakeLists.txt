set(HEADER_FILES
    optimizer.h
    types.h
)

# Install headers
install(FILES ${HEADER_FILES}
        DESTINATION include/OptimLight)