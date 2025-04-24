#include "Output.h"

Output::Output(Run *run)
{
    m_run = run;
}

void Output::Write2HDF5(const std::string &filename)
{
    auto N1local = static_cast<hsize_t>(m_run->m_mesh->getnumberCellsX());
    auto N2local = static_cast<hsize_t>(m_run->m_mesh->getnumberCellsY());
    auto N1global = static_cast<hsize_t>(m_run->m_mesh->getnumberCellsXGlobal());
    auto N2global = static_cast<hsize_t>(m_run->m_mesh->getnumberCellsYGlobal());
    auto offsetX = static_cast<hsize_t>(m_run->m_mesh->getoffsetX());
    auto offsetY = static_cast<hsize_t>(m_run->m_mesh->getoffsetY());

    H5::FileAccPropList fapl;
    H5Pset_fapl_mpio(fapl.getId(), MPI_COMM_WORLD, MPI_INFO_NULL);

    H5::H5File file(filename, H5F_ACC_TRUNC, H5::FileCreatPropList::DEFAULT, fapl);
    hsize_t global_dims[2] = {N1global, N2global};
    hsize_t local_dims[2] = {N1local, N2local};
    hsize_t offset[2] = {offsetX, offsetY};

    H5::DataSpace global_dataspace(2, global_dims);
    H5::DataSpace local_dataspace(2, local_dims);
    global_dataspace.selectHyperslab(H5S_SELECT_SET, local_dims, offset);

    H5::DSetCreatPropList plist;
    plist.setChunk(2, local_dims);

    H5::DataSet dataset_U = file.createDataSet("U", H5::PredType::NATIVE_DOUBLE, global_dataspace, plist);
    H5::DataSet dataset_V = file.createDataSet("V", H5::PredType::NATIVE_DOUBLE, global_dataspace, plist);
    H5::DataSet dataset_W = file.createDataSet("W", H5::PredType::NATIVE_DOUBLE, global_dataspace, plist);
    H5::DataSet dataset_P = file.createDataSet("P", H5::PredType::NATIVE_DOUBLE, global_dataspace, plist);
    H5::DataSet dataset_T = file.createDataSet("T", H5::PredType::NATIVE_DOUBLE, global_dataspace, plist);
    H5::DataSet dataset_Rho = file.createDataSet("Rho", H5::PredType::NATIVE_DOUBLE, global_dataspace, plist);
    H5::DataSet dataset_numchild = file.createDataSet("Numchild", H5::PredType::NATIVE_INT, global_dataspace, plist);
    H5::DataSet dataset_mfp = file.createDataSet("MFP", H5::PredType::NATIVE_DOUBLE, global_dataspace, plist);

    std::vector<double> U(N1local * N2local, 0.0);
    std::vector<double> V(N1local * N2local, 0.0);
    std::vector<double> W(N1local * N2local, 0.0);
    std::vector<double> P(N1local * N2local, 0.0);
    std::vector<double> T(N1local * N2local, 0.0);
    std::vector<double> Rho(N1local * N2local, 0.0);
    std::vector<int> Numchild(N1local * N2local, 0);
    std::vector<double> MFP (N1local * N2local, 0.0);
    for (size_t i = 0; i < N1local; ++i){
        for (size_t j = 0; j < N2local; ++j){
            auto cell = m_run->m_cells[i * N2local + j];
            auto phase = cell.getphase();
            auto element = cell.getelement();

            U[i * N2local + j] = phase->getvelocity()[0];
            V[i * N2local + j] = phase->getvelocity()[1];
            W[i * N2local + j] = phase->getvelocity()[2];
            P[i * N2local + j] = phase->getpressure();
            T[i * N2local + j] = phase->gettemperature();
            Rho[i * N2local + j] = phase->getdensity();
            Numchild[i * N2local + j] = cell.getchildren().size();
            // for(auto& child : cell.getchildren()){
            //     Numchild[i * N2local + j] += child->getchildren().size();
            // }
            MFP[i * N2local + j] = cell.getmfp();
        }
    }

    dataset_U.write(U.data(), H5::PredType::NATIVE_DOUBLE, local_dataspace, global_dataspace);
    dataset_V.write(V.data(), H5::PredType::NATIVE_DOUBLE, local_dataspace, global_dataspace);
    dataset_W.write(W.data(), H5::PredType::NATIVE_DOUBLE, local_dataspace, global_dataspace);
    dataset_P.write(P.data(), H5::PredType::NATIVE_DOUBLE, local_dataspace, global_dataspace);
    dataset_T.write(T.data(), H5::PredType::NATIVE_DOUBLE, local_dataspace, global_dataspace);
    dataset_Rho.write(Rho.data(), H5::PredType::NATIVE_DOUBLE, local_dataspace, global_dataspace);
    dataset_numchild.write(Numchild.data(), H5::PredType::NATIVE_INT, local_dataspace, global_dataspace);
    dataset_mfp.write(MFP.data(), H5::PredType::NATIVE_DOUBLE, local_dataspace, global_dataspace);
}

void Output::Write2VTK(const std::string &filename)
{
    auto N1local = m_run->m_mesh->getnumberCellsX();
    auto N2local = m_run->m_mesh->getnumberCellsY();
    auto N1global = m_run->m_mesh->getnumberCellsXGlobal();
    auto N2global = m_run->m_mesh->getnumberCellsYGlobal();
    auto offsetX = m_run->m_mesh->getoffsetX();
    auto offsetY = m_run->m_mesh->getoffsetY();
    auto dx = m_run->m_mesh->getUnidX();
    auto dy = m_run->m_mesh->getUnidX();
    auto myid = m_run->myid;
    auto numprocs = m_run->numprocs;

    auto& local_elements = m_run->m_mesh->getElements();
    // create a mesh obj
    vtkSmartPointer<vtkStructuredGrid> local_grid = vtkSmartPointer<vtkStructuredGrid>::New();
    local_grid->SetDimensions(N1local, N2local, 1);

    // create a points set
    vtkSmartPointer<vtkPoints> local_points = vtkSmartPointer<vtkPoints>::New();
    for(auto& element : local_elements){
        auto& vertices = element->getvertices();
        
        auto x = vertices[0]->getPosition()[0];
        auto y = vertices[0]->getPosition()[1];
        double z = 0.0;
        local_points->InsertNextPoint(x, y, z);
    }

    local_grid->SetPoints(local_points);

    std::vector<double> U(N1local * N2local, 0.0);
    std::vector<double> V(N1local * N2local, 0.0);
    std::vector<double> W(N1local * N2local, 0.0);
    std::vector<double> P(N1local * N2local, 0.0);
    std::vector<double> T(N1local * N2local, 0.0);
    std::vector<double> Rho(N1local * N2local, 0.0);
    for (size_t i = 0; i < N1local; ++i){
        for (size_t j = 0; j < N2local; ++j){
            auto cell = m_run->m_cells[i * N2local + j];
            auto phase = cell.getphase();
            auto element = cell.getelement();

            U[i * N2local + j] = phase->getvelocity()[0];
            V[i * N2local + j] = phase->getvelocity()[1];
            W[i * N2local + j] = phase->getvelocity()[2];
            P[i * N2local + j] = phase->getpressure();
            T[i * N2local + j] = phase->gettemperature();
            Rho[i * N2local + j] = phase->getdensity();
        }
    }

    // 创建并添加标量场数据
    auto addScalarField = [&](const std::vector<double>& data, const std::string& name) {
        vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
        array->SetName(name.c_str());
        array->SetNumberOfComponents(1);
        array->SetNumberOfTuples(N1local * N2local);
        for (size_t idx = 0; idx < data.size(); ++idx) {
            array->SetValue(idx, data[idx]);
        }
        local_grid->GetPointData()->AddArray(array);
    };

    
    addScalarField(U, "U");
    addScalarField(V, "V");
    addScalarField(P, "P");
    addScalarField(T, "T");
    addScalarField(Rho, "Rho");

    // // wtite in vtk
    // vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    // writer->SetFileName(filename.c_str());
    // // writer->SetDataModeToAscii(); 
    // writer->SetDataModeToBinary();  
    // writer->SetInputData(local_grid);
    // writer->Write();

    // set local grid scale
    int local_extent[6] = {
        offsetX, offsetX + static_cast<int>(N1local) - 1,
        offsetY, offsetY + static_cast<int>(N2local) - 1,
        0, 0
    };

    local_grid->SetExtent(local_extent);

   // Write the local .vts file
   vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
   writer->SetFileName((filename + "_" + std::to_string(myid) + ".vts").c_str());
   writer->SetDataModeToBinary();
   writer->SetInputData(local_grid);
   writer->Write();

   // Write the parallel .pvts file (only on the root process)
   if (myid == 0) {
        auto pwriter =  vtkSmartPointer<vtkXMLPStructuredGridWriter>::New();
        pwriter->SetFileName((filename + ".pvts").c_str());
        pwriter->SetNumberOfPieces(numprocs);
        pwriter->SetInputData(local_grid);
        pwriter->SetDataModeToBinary();
        pwriter->Update();
    }
}
