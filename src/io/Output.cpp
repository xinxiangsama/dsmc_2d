#include "Output.h"
#include <vtkUnstructuredGrid.h>
#include <vtkQuad.h>
#include <vtkCellData.h>
#include <vtkXMLUnstructuredGridWriter.h>

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
    H5::DataSet dataset_numgrandchild = file.createDataSet("Numgrandchild", H5::PredType::NATIVE_INT, global_dataspace, plist);
    H5::DataSet dataset_mfp = file.createDataSet("MFP", H5::PredType::NATIVE_DOUBLE, global_dataspace, plist);

    std::vector<double> U(N1local * N2local, 0.0);
    std::vector<double> V(N1local * N2local, 0.0);
    std::vector<double> W(N1local * N2local, 0.0);
    std::vector<double> P(N1local * N2local, 0.0);
    std::vector<double> T(N1local * N2local, 0.0);
    std::vector<double> Rho(N1local * N2local, 0.0);
    std::vector<int> Numchild(N1local * N2local, 0);
    std::vector<int> Numgrandchild(N1local * N2local, 0);
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
            for(auto& child : cell.getchildren()){
                Numgrandchild[i * N2local + j] += child->getchildren().size();
            }
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
    dataset_numgrandchild.write(Numgrandchild.data(), H5::PredType::NATIVE_INT, local_dataspace, global_dataspace);
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
    // writer->SetDataModeToBinary();
    writer->SetDataModeToAscii(); 
    writer->SetInputData(local_grid);
    writer->Write();

   // Write the parallel .pvts file (only on the root process)
   if (myid == 0) {
        WriteParallelVTSHeader(filename, numprocs, N1global, N2global, N1local, N2local);
    }
}


void Output::flattenCells(Cell* cell, std::vector<Cell*>& flatlist)
{
    flatlist.push_back(cell);
    for (auto& child : cell->getchildren()) {
        flattenCells(child.get(), flatlist); // 递归添加子单元格
    }
}

void Output::WriteAMRmesh(const std::string &filename)
{
    std::vector<Cell*> all_cells;
    for (auto& root_cell : m_run->m_cells) {
        flattenCells(&root_cell, all_cells); // Flatten from each top-level cell
    }

    hsize_t total_cells = static_cast<hsize_t>(all_cells.size());

    H5::FileAccPropList fapl;
    H5Pset_fapl_mpio(fapl.getId(), MPI_COMM_WORLD, MPI_INFO_NULL);
    H5::H5File file(filename, H5F_ACC_TRUNC, H5::FileCreatPropList::DEFAULT, fapl);

    hsize_t dims[1] = {total_cells};
    H5::DataSpace dataspace(1, dims);

    std::vector<double> X(total_cells), Y(total_cells), Dx(total_cells), Dy(total_cells), Mfp(total_cells);
    std::vector<int> Level(total_cells);

    for (size_t i = 0; i < total_cells; ++i) {
        Cell* cell = all_cells[i];
        auto elem = cell->getelement();
        auto pos = elem->getposition();
        X[i] = pos[0];
        Y[i] = pos[1];
        Dx[i] = elem->getL1();
        Dy[i] = elem->getL2();
        Mfp[i] = cell->getmfp();
        Level[i] = static_cast<int>(cell->getAMRlevel()); // Lv1 = 0, Lv2 = 1, ...
    }

    file.createDataSet("X", H5::PredType::NATIVE_DOUBLE, dataspace).write(X.data(), H5::PredType::NATIVE_DOUBLE);
    file.createDataSet("Y", H5::PredType::NATIVE_DOUBLE, dataspace).write(Y.data(), H5::PredType::NATIVE_DOUBLE);
    file.createDataSet("Dx", H5::PredType::NATIVE_DOUBLE, dataspace).write(Dx.data(), H5::PredType::NATIVE_DOUBLE);
    file.createDataSet("Dy", H5::PredType::NATIVE_DOUBLE, dataspace).write(Dy.data(), H5::PredType::NATIVE_DOUBLE);
    file.createDataSet("MFP", H5::PredType::NATIVE_DOUBLE, dataspace).write(Mfp.data(), H5::PredType::NATIVE_DOUBLE);
    file.createDataSet("Level", H5::PredType::NATIVE_INT, dataspace).write(Level.data(), H5::PredType::NATIVE_INT);
}

void Output::WriteAMR2VTK(const std::string &filename)
{
    std::vector<Cell*> all_cells;
    for (auto& root_cell : m_run->m_cells) {
        flattenCells(&root_cell, all_cells); // Flatten AMR tree
    }

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // 存储 Cell 的物理量
    auto MFPArray = vtkSmartPointer<vtkDoubleArray>::New();
    MFPArray->SetName("MFP");
    MFPArray->SetNumberOfComponents(1);

    auto LevelArray = vtkSmartPointer<vtkIntArray>::New();
    LevelArray->SetName("Level");
    LevelArray->SetNumberOfComponents(1);

    vtkIdType pointId = 0;
    for (auto* cell : all_cells)
    {
        auto elem = cell->getelement();
        auto pos = elem->getposition();
        double cx = pos[0], cy = pos[1];
        double dx = elem->getL1(), dy = elem->getL2();

        // 计算四个角点（逆时针）
        std::array<std::array<double, 3>, 4> corners {{
            {cx - dx/2, cy - dy/2, 0.0},
            {cx + dx/2, cy - dy/2, 0.0},
            {cx + dx/2, cy + dy/2, 0.0},
            {cx - dx/2, cy + dy/2, 0.0}
        }};

        std::array<vtkIdType, 4> pt_ids;
        for (int i = 0; i < 4; ++i) {
            pt_ids[i] = points->InsertNextPoint(corners[i].data());
        }

        // 构建 Quad 单元格
        auto quad = vtkSmartPointer<vtkQuad>::New();
        for (int i = 0; i < 4; ++i) {
            quad->GetPointIds()->SetId(i, pt_ids[i]);
        }
        grid->InsertNextCell(quad->GetCellType(), quad->GetPointIds());

        // 添加属性
        MFPArray->InsertNextValue(cell->getmfp());
        LevelArray->InsertNextValue(static_cast<int>(cell->getAMRlevel()));
    }

    // 设置点和属性
    grid->SetPoints(points);
    grid->GetCellData()->AddArray(MFPArray);
    grid->GetCellData()->AddArray(LevelArray);

    // 写入 .vtu 文件
    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName((filename + ".vtu").c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToBinary();
    // writer->SetDataModeToAscii();
    writer->Write();
}

void Output::WriteParallelVTSHeader(const std::string &filename, int numprocs, int N1global, int N2global, int N1local, int N2local)
{
    std::ofstream file(filename + ".pvts");
    file << R"(<?xml version="1.0"?>)" << "\n";
    file << R"(<VTKFile type="PStructuredGrid" version="0.1" byte_order="LittleEndian">)" << "\n";
    file << "  <PStructuredGrid WholeExtent=\"0 " << (N1global - 1)
         << " 0 " << (N2global - 1) << " 0 0\" GhostLevel=\"0\">" << "\n";
    file << "    <PPoints>" << "\n";
    file << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\"/>" << "\n";
    file << "    </PPoints>" << "\n";

    file << "    <PPointData>" << "\n";
    file << "      <PDataArray type=\"Float64\" Name=\"U\"/>" << "\n";
    file << "      <PDataArray type=\"Float64\" Name=\"V\"/>" << "\n";
    file << "      <PDataArray type=\"Float64\" Name=\"P\"/>" << "\n";
    file << "      <PDataArray type=\"Float64\" Name=\"T\"/>" << "\n";
    file << "      <PDataArray type=\"Float64\" Name=\"Rho\"/>" << "\n";
    file << "    </PPointData>" << "\n";

    for (int pid = 0; pid < numprocs; ++pid) {
        int i_start = (pid % (N1global / N1local)) * N1local;
        int j_start = (pid / (N1global / N1local)) * N2local;

        int i_end = i_start + N1local - 1;
        int j_end = j_start + N2local - 1;

        file << "    <Piece Extent=\"" << i_start << " " << i_end << " "
             << j_start << " " << j_end << " 0 0\" Source=\""
             << filename << "_" << pid << ".vts\"/>" << "\n";
    }

    file << "  </PStructuredGrid>" << "\n";
    file << "</VTKFile>" << std::endl;
}
